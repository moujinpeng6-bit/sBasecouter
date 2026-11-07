"""Convert verified spatial assets into one AnnData record per slice and bin size."""

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import anndata
import numpy as np

from . import __version__
from .formats import expression, spatial

SUPPORTED = {"Visium", "Visium HD", "Xenium"}
HUMAN = {"homo sapiens", "human"}


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def checked_files(row, report_path):
    assets = row["files"]
    needed = {"image", "cells", "experiment"} if row["platform"] == "Xenium" else {"image", "positions", "scalefactors"}
    if "h5" in assets or row.get("status") == "recount_ready_for_validation":
        needed.add("h5")
    else:
        needed.update(("matrix", "features", "barcodes"))
    paths, provenance = {}, {}
    for role in sorted(needed):
        asset = assets.get(role)
        if role == "image" and not asset and row["platform"] != "Xenium":
            asset = assets.get("image_hires")
        if not asset:
            raise ValueError(f"Report has no {role} asset")
        path = (report_path.parent / asset["path"]).resolve(strict=True)
        digest = sha256(path)
        if digest != asset.get("sha256"):
            raise ValueError(f"{role} SHA-256 differs from verified report")
        paths[role] = path
        provenance[role] = {"path": str(path), "sha256": digest, "bytes": path.stat().st_size}
        provenance[role].update({key:asset[key] for key in ("source_url","member") if key in asset})
    return paths, provenance


def write_record(row, request, report_path, report_sha, destination, normalize_total):
    paths, inputs = checked_files(row, report_path)
    counts, barcodes, features = expression(paths)
    obs, native_xy, image_xy, frame = spatial(paths, row["platform"], barcodes, row.get("bin_um"))
    obs["total_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    obs["detected_features"] = np.diff(counts.indptr)
    gene_mask = features["feature_type"].eq("Gene Expression").to_numpy()
    gene_total = (obs["total_counts"].to_numpy() if gene_mask.all()
                  else np.asarray(counts[:, gene_mask].sum(axis=1)).ravel())
    obs["gene_expression_counts"] = gene_total
    data = anndata.AnnData(X=counts, obs=obs, var=features)
    data.obsm["spatial"] = native_xy
    data.obsm["image_xy"] = image_xy
    if normalize_total is not None:
        if not gene_mask.any():
            raise ValueError("No Gene Expression features for optional normalization")
        factors = np.divide(normalize_total, gene_total,
                            out=np.zeros(len(obs), dtype=np.float64), where=gene_total > 0)
        log_counts = counts.astype(np.float64).multiply(gene_mask).multiply(factors[:, None]).tocsr()
        log_counts.eliminate_zeros()
        np.log1p(log_counts.data, out=log_counts.data)
        data.layers["log1p_total"] = log_counts.astype(np.float32)

    sample = {key:row[key] for key in ("sample_id","source_sample_id","study_id","platform","source","source_url",
              "condition_evidence","clinical_status","observed_sample_state","lineage") if key in row}
    if not sample.get("observed_sample_state"):
        sample["observed_sample_state"] = (row.get("condition_evidence") or {}).get("sample_state","unknown")
    provenance = {"sample":sample,"search_request":request,"input_files":inputs,
                  "report_sha256":report_sha,"matrix_kind":"publisher_processed_counts",
                  "builder_version":__version__}
    if row["status"] == "recount_ready_for_validation":
        provenance.update(matrix_kind="local_recount", recount=row["recount"])
    frame["slice_id"] = row["sample_id"]
    data.uns["spatial_reference"] = frame
    data.uns["source_record_json"] = json.dumps(provenance, ensure_ascii=False)
    if normalize_total is not None:
        data.uns["normalization"] = {"layer":"log1p_total","target_sum":normalize_total,
            "denominator":"Gene Expression counts per observation","non_gene_features":"zero in layer, preserved in X"}

    with tempfile.NamedTemporaryFile(dir=destination.parent, suffix=".h5ad", delete=False) as stream:
        temporary = Path(stream.name)
    try:
        data.write_h5ad(temporary, compression="gzip")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return {"file":destination.name,"sha256":sha256(destination),"observations":data.n_obs,"features":data.n_vars,
            "gene_expression_features":int(gene_mask.sum()),"nonzero_counts":counts.nnz,
            "spatial_reference":frame,"source_record":provenance}


def eligibility(row, request, *, local_recount=False):
    required = "recount_ready_for_validation" if local_recount else "correspondence_passed"
    if row.get("status") != required:
        return "source_spatial_validation_not_passed"
    if row.get("platform") not in SUPPORTED:
        return "platform_not_supported"
    if row.get("matrix_scope") == "header_only":
        return "full_count_matrix_required"
    declared = [str(value).casefold().strip() for value in (request.get("organism"), row.get("organism")) if value and value!="*"]
    if not declared or any(value not in HUMAN for value in declared):
        return "human_species_not_confirmed"
    return None


def build(reports, out, *, sample_ids=(), normalize_total=None, local_recount=False):
    if normalize_total is not None and (not math.isfinite(normalize_total) or normalize_total <= 0):
        raise ValueError("Normalization target must be positive and finite")
    wanted = set(sample_ids)
    samples = {}
    for path in reports:
        path = path.resolve(strict=True)
        body = path.read_bytes()
        report = json.loads(body)
        request = report.get("request", {})
        report_sha = hashlib.sha256(body).hexdigest()
        for row in report["samples"]:
            sample = row["sample_id"]
            if wanted and sample not in wanted:
                continue
            if sample in samples:
                raise ValueError(f"Duplicate sample ID across reports: {sample}")
            samples[sample] = (row, request, path, report_sha)
    if wanted - samples.keys():
        raise ValueError("Unknown sample IDs: " + ", ".join(sorted(wanted - samples.keys())))
    out.mkdir(parents=True, exist_ok=False)
    catalog = {"schema":1,"builder_version":__version__,"normalization_target":normalize_total,"samples":[]}

    def save():
        pending = out / "catalog.tmp"
        with pending.open("w", encoding="utf8") as stream:
            json.dump(catalog, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        pending.replace(out / "catalog.json")

    save()
    for sample, (row, request, path, report_sha) in samples.items():
        result = {key:row[key] for key in ("sample_id","source_sample_id","study_id","platform","source_url") if key in row}
        reason = eligibility(row, request, local_recount=local_recount)
        if reason:
            result.update(status="skipped", reason=reason)
        else:
            filename = hashlib.sha256(sample.encode()).hexdigest() + ".h5ad"
            try:
                result.update(write_record(row, request, path, report_sha, out/filename, normalize_total))
                result["status"] = "built"
            except (OSError, EOFError, ValueError, KeyError, TypeError, IndexError, OverflowError) as exc:
                result.update(status="failed", reason=str(exc), error_type=type(exc).__name__)
        catalog["samples"].append(result)
        save()
    catalog["totals"] = {status:sum(row["status"]==status for row in catalog["samples"])
                         for status in ("built","failed","skipped")}
    save()
    return catalog
