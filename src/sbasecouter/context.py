"""Create a spatial expression layer without changing observed counts."""

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import anndata
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from .build import sha256

UNKNOWN = "unknown"


def neighbor_graph(xy, neighbors, radius):
    n = len(xy)
    k = min(neighbors, n - 1)
    if k <= 0:
        return sparse.csr_matrix((n, n), dtype=np.float32)
    tree = cKDTree(xy)
    pieces = []
    for start in range(0, n, 16384):
        stop = min(n, start + 16384)
        distances, indices = tree.query(xy[start:stop], k=k + 1,
                                        distance_upper_bound=np.nextafter(radius, np.inf))
        valid = (indices != np.arange(start, stop)[:, None]) & (indices < n) & (distances <= radius)
        # Duplicate coordinates can cause the queried cell to be absent even with k+1 candidates.
        valid &= np.cumsum(valid, axis=1) <= k
        rows, cols = np.where(valid)
        pieces.append(sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, indices[rows, cols])),
                                        shape=(stop - start, n)))
    return sparse.vstack(pieces, format="csr")


def cell_labels(data, label_file, slice_id):
    if label_file is None:
        if "cell_type" in data.obs:
            return data.obs["cell_type"].astype("string").fillna(UNKNOWN).str.strip().replace("", UNKNOWN), "obs.cell_type"
        return pd.Series(UNKNOWN, index=data.obs_names, dtype="string"), "unavailable"
    table = pd.read_csv(label_file, dtype=str, keep_default_na=False)
    if not {"slice_id", "barcode", "cell_type"}.issubset(table.columns):
        raise ValueError("Cell type CSV requires slice_id, barcode and cell_type columns")
    if not table["slice_id"].eq(slice_id).all():
        raise ValueError("Cell type CSV slice ID differs from the input slice")
    if table["barcode"].eq("").any() or table["barcode"].duplicated().any():
        raise ValueError("Cell type CSV contains an empty or duplicate barcode")
    extra = pd.Index(table["barcode"]).difference(data.obs_names)
    if len(extra) or not len(table):
        raise ValueError("Cell type CSV must contain matching barcodes from this slice only")
    labels = table.set_index("barcode")["cell_type"].reindex(data.obs_names)
    return labels.astype("string").fillna(UNKNOWN).str.strip().replace("", UNKNOWN), "matched_cell_type_csv"


def neighborhood_coordinates(data, frame):
    unit = frame.get("observation_unit")
    if unit == "cell" and frame.get("native_unit") == "micrometre":
        return np.asarray(data.obsm["spatial"], dtype=np.float64), np.arange(data.n_obs), 50.0, "cell_centroids"
    fields = ["array_col", "array_row", "in_tissue"]
    if unit not in {"spot", "bin"} or not set(fields).issubset(data.obs):
        raise ValueError("Expected cell centroids in micrometres or validated spot/bin array indices")
    grid = data.obs[fields[:2]].to_numpy(dtype=np.float64)
    tissue = data.obs["in_tissue"].to_numpy()
    if (not np.isfinite(grid).all() or (grid < 0).any() or (grid != np.floor(grid)).any()
            or not np.isin(tissue, [0, 1]).all() or data.obs.duplicated(fields[:2]).any()):
        raise ValueError("Expected unique, nonnegative array positions and binary in_tissue flags")
    active = np.flatnonzero(tissue)
    if unit == "spot":
        if np.any((grid[:, 0] - grid[:, 1]) % 2):
            raise ValueError("Visium spot rows and columns disagree with the staggered layout")
        return grid * [50, 50 * math.sqrt(3)], active, 150.0, "visium_staggered_capture_grid"
    size = float(frame["bin_um"])
    if not math.isfinite(size) or size <= 0:
        raise ValueError("Expected a positive native Visium HD bin size")
    return grid * size, active, 2 * size, "visium_hd_square_capture_grid"


def contextualize(data, *, neighbors, radius, alpha, same_type_weight, label_file=None):
    frame = data.uns.get("spatial_reference", {})
    if not isinstance(frame.get("slice_id"), str) or not frame["slice_id"]:
        raise ValueError("Expected a named, validated slice")
    unit = frame.get("observation_unit")
    if unit != "cell" and label_file is not None:
        raise ValueError("Cell type labels require cell-level observations")
    if "slice_id" in data.obs and not data.obs["slice_id"].eq(frame["slice_id"]).all():
        raise ValueError("A context run accepts one slice only")
    if "spatial_context" in data.uns:
        raise ValueError("Use the original standardized h5ad to compute a new context")
    if not data.n_obs or not data.obs_names.is_unique or (data.obs_names == "").any():
        raise ValueError("Expected unique observation barcodes")
    xy, active, default_radius, geometry = neighborhood_coordinates(data, frame)
    radius = default_radius if radius is None else radius
    if xy.shape != (data.n_obs, 2) or not np.isfinite(xy).all():
        raise ValueError("Expected finite 2D neighborhood coordinates")
    if not 1 <= neighbors <= 128 or not math.isfinite(radius) or radius <= 0:
        raise ValueError("Expected 1–128 neighbors and a positive finite radius")
    if not math.isfinite(alpha) or not 0 <= alpha <= 1 or not math.isfinite(same_type_weight) or not 1 <= same_type_weight <= 10:
        raise ValueError("Expected alpha in [0, 1] and a same-type weight in [1, 10]")

    counts = sparse.csr_matrix(data.X)
    values = counts.data
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Expected nonnegative integer observed counts in X")
    genes = data.var["feature_type"].eq("Gene Expression").to_numpy(dtype=np.float32)
    if not genes.any():
        raise ValueError("No Gene Expression features in the input")
    rna = counts.astype(np.float32).multiply(genes).tocsr()
    rna.eliminate_zeros()
    totals = np.asarray(rna.sum(axis=1)).ravel()
    has_rna = totals > 0
    rates = rna.multiply(np.divide(10000, totals, out=np.zeros_like(totals), where=has_rna)[:, None]).tocsr()

    graph = neighbor_graph(xy[active], neighbors, radius)
    if len(active) != data.n_obs:
        local = graph.tocoo()
        graph = sparse.csr_matrix((local.data, (active[local.row], active[local.col])),
                                  shape=(data.n_obs, data.n_obs))
    degree = np.diff(graph.indptr)
    cols = graph.indices
    weights = has_rna[cols].astype(np.float32)
    label_source, categories, known_count = "not_applicable", [], 0
    if unit == "cell":
        labels, label_source = cell_labels(data, label_file, frame["slice_id"])
        labels = labels.mask(labels.str.casefold().eq(UNKNOWN), UNKNOWN)
        codes, categories = pd.factorize(labels, sort=True)
        unknown_code = categories.get_loc(UNKNOWN) if UNKNOWN in categories else -1
        known = codes != unknown_code
        known_count = int(known.sum())
        rows = np.repeat(np.arange(data.n_obs), degree)
        membership = sparse.csr_matrix((np.ones(data.n_obs, dtype=np.float32), (np.arange(data.n_obs), codes)),
                                       shape=(data.n_obs, len(categories)))
        inverse = np.divide(1, degree, out=np.zeros(data.n_obs, dtype=np.float32), where=degree > 0)
        data.obsm["neighbor_type_fractions"] = (graph @ membership).multiply(inverse[:, None]).tocsr()
        weights[known[rows] & (codes[rows] == codes[cols])] *= same_type_weight
        data.obs["cell_type"] = pd.Categorical(labels)
    eligible = sparse.csr_matrix((weights, graph.indices, graph.indptr), shape=graph.shape, copy=True)
    weight_sum = np.asarray(eligible.sum(axis=1)).ravel()
    available = weight_sum > 0
    eligible.eliminate_zeros()
    rna_neighbor_count = np.diff(eligible.indptr)
    normalized = eligible.multiply(np.divide(1, weight_sum, out=np.zeros_like(weight_sum),
                                              where=available)[:, None]).tocsr()
    neighborhood = normalized @ rates
    inferred = ~has_rna & available & (alpha > 0)
    blend = np.where(inferred, 1, np.where(available, alpha, 0)).astype(np.float32)
    context = (rates.multiply((1 - blend)[:, None]) + neighborhood.multiply(blend[:, None])).tocsr()

    parameters = {"k":neighbors, "radius_um":radius, "alpha":alpha, "same_type_weight":same_type_weight}
    data.obs["spatial_observation_id"] = frame["slice_id"] + "/" + data.obs_names
    if unit == "cell":
        data.obs["spatial_cell_id"] = data.obs["spatial_observation_id"]
    data.obs["spatial_neighbors"] = degree
    data.obs["rna_neighbors"] = rna_neighbor_count
    data.obs["context_imputed_from_neighbors_only"] = inferred
    data.obsp["spatial_neighbors"] = graph
    data.obsp["spatial_rna_weights"] = normalized
    data.obsm["context_xy_um"] = xy
    data.layers["neighbor_expression"] = neighborhood
    data.layers["spatial_context_expression"] = context
    data.uns["spatial_context"] = {**parameters, "cell_type_source":label_source, "neighborhood_geometry":geometry,
        "neighbor_type_columns":list(categories), "expression_scale":"RNA counts per 10000, linear",
        "method":"radius-limited kNN, label-weighted neighbor mean and explicit self-neighbor blend",
        "expression_semantics":"derived estimate; X contains the original observed counts"}
    counts_summary = {"observations":data.n_obs, "observations_with_neighbors":int((degree > 0).sum()),
                      "observations_with_rna_neighbors":int(available.sum())}
    cell_summary = {"cells":data.n_obs, "cells_with_neighbors":counts_summary["observations_with_neighbors"],
                    "cells_with_rna_neighbors":counts_summary["observations_with_rna_neighbors"],
                    "cells_with_known_types":known_count} if unit == "cell" else {}
    return {"slice_id":frame["slice_id"], "observation_unit":unit, **counts_summary, **cell_summary,
            "rna_genes":int(genes.sum()), "cell_type_source":label_source,
            "neighbor_only_estimates":int(inferred.sum()), "parameters":parameters}


def build_context(source, destination, *, neighbors, radius, alpha, same_type_weight, label_file=None, label_source=None):
    source, destination = source.resolve(strict=True), destination.resolve()
    if destination.suffix.lower() != ".h5ad" or destination.exists():
        raise ValueError("Choose a new .h5ad output file")
    if label_source and label_file is None:
        raise ValueError("A cell type source requires a cell type CSV")
    destination.parent.mkdir(parents=True, exist_ok=True)
    data = anndata.read_h5ad(source)
    result = contextualize(data, neighbors=neighbors, radius=radius, alpha=alpha,
                           same_type_weight=same_type_weight, label_file=label_file)
    provenance = data.uns["spatial_context"]
    provenance["input_sha256"] = sha256(source)
    if label_file is not None:
        provenance["labels_sha256"] = sha256(label_file)
        if label_source:
            provenance["labels_reference"] = label_source
    run_id = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    provenance["run_id"] = run_id
    data.obs["spatial_context_id"] = data.obs["spatial_observation_id"] + "/" + run_id
    with tempfile.NamedTemporaryFile(dir=destination.parent, suffix=".h5ad", delete=False) as file:
        temporary = Path(file.name)
    try:
        data.write_h5ad(temporary, compression="gzip")
        # Reserve the output before replacing it, so a concurrent run cannot overwrite a result.
        fd = os.open(destination, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
        try:
            temporary.replace(destination)
        except OSError:
            destination.unlink()
            raise
    finally:
        temporary.unlink(missing_ok=True)
    return {**result, "output":str(destination)}
