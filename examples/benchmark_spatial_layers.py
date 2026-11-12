"""Compare expression features against independent layer labels by held-out donor."""

import argparse
from hashlib import file_digest
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from sbasecouter.context import neighbor_graph, neighborhood_coordinates
from sbasecouter.formats import expression, spatial


def read_slice(root, row):
    paths = {role: root / name for role, name in row["files"].items()}
    for role, path in {**paths, "labels": root / row["labels"]}.items():
        with path.open("rb") as stream:
            if file_digest(stream, "sha256").hexdigest() != row["sha256"][role]:
                raise ValueError(f"Changed {role} file for {row['slice']}")
    counts, barcodes, features = expression(paths)
    obs, _, _, frame = spatial(paths, "Visium", barcodes)
    annotation = pd.read_csv(root / row["labels"], dtype=str)
    if annotation["barcode"].duplicated().any():
        raise ValueError(f"Duplicate annotation barcode for {row['slice']}")
    labels = annotation.set_index("barcode")["layer"].reindex(barcodes).fillna("").to_numpy()
    counts = counts[:, features["feature_type"].eq("Gene Expression").to_numpy()].astype(np.float32)
    features = features.loc[features["feature_type"].eq("Gene Expression"), "feature_id"]
    depth = np.asarray(counts.sum(axis=1)).ravel()
    observed = depth > 0
    rates = counts.multiply(np.divide(10000, depth, out=np.zeros_like(depth), where=observed)[:, None]).tocsr()
    xy, active, radius, _ = neighborhood_coordinates(SimpleNamespace(obs=obs), frame)
    tissue = active[observed[active]]
    graph = neighbor_graph(xy[tissue], 8, radius)
    degree = np.diff(graph.indptr)
    weights = graph.multiply(np.divide(1, degree, out=np.zeros_like(degree, dtype=np.float32), where=degree > 0)[:, None]).tocsr()
    scored = labels[tissue] != ""
    if not scored.any() or not features.is_unique:
        raise ValueError(f"No usable labels or duplicate genes in {row['slice']}")
    row = {**row, "annotations_supplied": len(annotation), "annotations_scored": int(scored.sum())}
    return row, pd.Index(features), rates[tissue], weights, degree, labels[tissue], scored


def feature_statistics(matrix):
    logged = matrix.copy()
    np.log1p(logged.data, out=logged.data)
    return np.asarray(logged.sum(axis=0)).ravel(), np.asarray(logged.power(2).sum(axis=0)).ravel()


def design(item, columns):
    _, _, rates, graph, degree, labels, scored = item
    current = rates[:, columns].toarray()
    mixed = graph @ current
    mixed = np.where((degree > 0)[:, None], 0.75 * current + 0.25 * mixed, current)
    return np.log1p(current[scored]), np.log1p(mixed[scored]), labels[scored]


def score(train, labels, test, truth):
    classes = np.unique(labels)
    if not np.isin(truth, classes).all():
        raise ValueError("Held-out donor has a label absent from training")
    means = np.stack([train[labels == label].mean(axis=0) for label in classes])
    precision = 1 / np.maximum(train.var(axis=0), 0.05)
    weighted = means * precision
    predicted = classes[(test @ weighted.T - (0.5 * means * weighted).sum(axis=1)).argmax(axis=1)]
    f1 = []
    for label in np.unique(truth):
        actual, guess = truth == label, predicted == label
        f1.append(2 * np.count_nonzero(actual & guess) / (actual.sum() + guess.sum()))
    return {"accuracy": float((predicted == truth).mean()), "macro_f1": float(np.mean(f1))}


def benchmark(manifest):
    specification = json.loads(manifest.read_text(encoding="utf8"))
    records = specification["samples"]
    if len({row["slice"] for row in records}) != len(records):
        raise ValueError("Slice IDs must be unique")
    cohorts = {row["donor"]: row["cohort"] for row in records}
    if any(cohorts[row["donor"]] != row["cohort"] for row in records):
        raise ValueError("Conflicting cohort labels for a donor")
    if set(cohorts.values()) != {"control", "schizophrenia"} or min(list(cohorts.values()).count(c) for c in set(cohorts.values())) < 2:
        raise ValueError("Both cohorts require at least two distinct donors")
    items = [read_slice(manifest.parent, row) for row in records]
    genes = items[0][1]
    for item in items[1:]:
        genes = genes.intersection(item[1], sort=False)
    if not len(genes):
        raise ValueError("No shared genes")
    items = [(row, genes, rates[:, index.get_indexer(genes)], graph, degree, labels, mask)
             for row, index, rates, graph, degree, labels, mask in items]
    statistics = [feature_statistics(item[2][item[6]]) for item in items]
    folds = []
    for donor, cohort in sorted(cohorts.items()):
        training = [i for i, item in enumerate(items) if item[0]["donor"] != donor]
        held_out = [i for i, item in enumerate(items) if item[0]["donor"] == donor]
        n = sum(items[i][6].sum() for i in training)
        mean = sum(statistics[i][0] for i in training) / n
        variance = sum(statistics[i][1] for i in training) / n - mean**2
        columns = np.argsort(-variance, kind="stable")[:512]
        train = [design(items[i], columns) for i in training]
        test = [design(items[i], columns) for i in held_out]
        train_labels, truth = [np.concatenate([x[2] for x in group]) for group in (train, test)]
        metrics = {name: score(np.concatenate([x[j] for x in train]), train_labels,
                               np.concatenate([x[j] for x in test]), truth)
                   for j, name in enumerate(("self", "with_neighbors"))}
        folds.append({"donor": donor, "cohort": cohort, "spots": len(truth), **metrics})
    summaries = {}
    for cohort in ("all", "control", "schizophrenia"):
        group = [fold for fold in folds if cohort == "all" or fold["cohort"] == cohort]
        summaries[cohort] = {name: {metric: float(np.mean([fold[name][metric] for fold in group]))
                                    for metric in ("accuracy", "macro_f1")}
                             for name in ("self", "with_neighbors")}
        summaries[cohort]["donors_with_improved_macro_f1"] = sum(
            fold["with_neighbors"]["macro_f1"] > fold["self"]["macro_f1"] for fold in group)
    return {"source": specification["source"], "label_source": specification["label_source"],
            "design": "leave one donor out; top 512 training-only variable genes; diagonal centroid classifier; mean across donors",
            "neighbors": {"maximum": 8, "radius_um": 150, "fraction": 0.25}, "shared_genes": len(genes),
            "slices": len(items), "folds": folds, "mean": summaries,
            "annotation_coverage": {item[0]["slice"]: {key: item[0][key] for key in
                                     ("annotations_supplied", "annotations_scored")} for item in items}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(benchmark(args.manifest), indent=2) + "\n", encoding="utf8")
