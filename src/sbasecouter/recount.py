"""Run vendor recounts and pass each spatial output through the existing builder."""

import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

from .build import HUMAN, build, sha256


RUN_ID = "spatial_count"
ASSAYS = {"Visium": {"three_prime", "probe"}, "Visium HD": {"three_prime", "probe"},
          "Xenium": {"in_situ"}}
POLICIES = {"three_prime": {"exonic": "--include-introns=false", "exonic_intronic": "--include-introns=true"},
            "probe": {"filtered_probes": "--filter-probes=true", "all_nondeprecated_probes": "--filter-probes=false"},
            "in_situ": {"resegment": None}}


def required_text(record, key):
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Manifest requires a nonempty {key}")
    return value.strip()


def asset(root, record, key, *, directory=False):
    path = (root / required_text(record, key)).resolve(strict=True)
    if not (path.is_dir() if directory else path.is_file()):
        raise ValueError(f"Invalid input {key}: {path}")
    return path


def fastq_pairs(directory, sample):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
        raise ValueError("Expected one Illumina FASTQ sample name")
    pattern = re.compile(rf"{re.escape(sample)}_(S\d+(?:_L\d{{3}})?)_([RI][12])_(\d{{3}})\.fastq\.gz")
    groups = {}
    inputs = []
    for path in sorted(directory.rglob("*.fastq.gz")):
        match = pattern.fullmatch(path.name)
        if not match or path.parent != directory or not path.is_file():
            raise ValueError("Use a flat FASTQ directory containing only the requested sample")
        lane, read, chunk = match.groups()
        groups.setdefault((lane, chunk), set()).add(read)
        inputs.append(path)
    if not groups or any(not {"R1", "R2"} <= reads for reads in groups.values()):
        raise ValueError("Expected complete R1/R2 FASTQ pairs for every lane and chunk")
    return inputs


def link_source(manifest, root, sample):
    if "source_report" not in manifest:
        return {"link": "manifest_declared"}
    path = asset(root, manifest, "source_report")
    body = path.read_bytes()
    report = json.loads(body)
    matches = [row for row in report["samples"] if row.get("sample_id") == sample["sample_id"]]
    if len(matches) != 1:
        raise ValueError("Expected exactly one matching sample in the SRAgent-spatial report")
    row = matches[0]
    for key in ("source_sample_id", "study_id", "platform"):
        if row.get(key) != sample[key]:
            raise ValueError(f"SRAgent-spatial report disagrees on {key}")
    declared = [value for value in (row.get("organism"), report.get("request", {}).get("organism"))
                if value and value != "*"]
    if any(value.casefold().strip() not in HUMAN for value in declared):
        raise ValueError("SRAgent-spatial report disagrees on organism")
    for key in ("source_url", "observed_sample_state", "condition_evidence", "lineage"):
        if key in row:
            sample[key] = row[key]
    return {"link": "sragent_spatial_report", "report_sha256": hashlib.sha256(body).hexdigest(),
            "sample_id": sample["sample_id"]}


def focus_image(bundle):
    names = ("morphology_focus/ch0000_dapi.ome.tif", "morphology_focus/morphology_focus_0000.ome.tif",
             "morphology_focus.ome.tif")
    return next((bundle / name for name in names if (bundle / name).is_file()), bundle / names[-1])


def create_plan(manifest, root, *, executable=None, cores=None, memory_gb=None):
    sample = {key: required_text(manifest, key) for key in
              ("sample_id", "source_sample_id", "study_id", "source_url", "organism", "platform")}
    if sample["organism"].casefold() not in HUMAN:
        raise ValueError("Only human samples are supported")
    for key in ("observed_sample_state", "condition_evidence", "lineage"):
        if key in manifest:
            sample[key] = manifest[key]
    platform = sample["platform"]
    assay = required_text(manifest, "assay")
    policy = required_text(manifest, "count_policy")
    if assay not in ASSAYS.get(platform, set()) or policy not in POLICIES[assay]:
        raise ValueError("Unsupported platform, assay or count policy combination")
    version = required_text(manifest, "engine_version")
    if not re.fullmatch(r"\d+\.\d+\.\d+", version) or tuple(map(int, version.split('.'))) < (4, 0, 0):
        raise ValueError("Pin a vendor engine version at 4.0.0 or later")
    if platform == "Xenium":
        engine = "xeniumranger"
        bundle = asset(root, manifest, "xenium_bundle", directory=True)
        for name in ("experiment.xenium", "morphology.ome.tif", "metrics_summary.csv", "analysis_summary.html"):
            if not (bundle / name).is_file():
                raise ValueError(f"Xenium resegmentation requires {name}")
        for name in ("transcripts", "cells"):
            if not ((bundle / f"{name}.zarr.zip").is_file() or (bundle / f"{name}.zarr").is_dir()):
                raise ValueError(f"Xenium resegmentation requires {name}.zarr.zip or {name}.zarr")
        if not focus_image(bundle).is_file():
            raise ValueError("Xenium resegmentation requires a morphology focus image")
        command = ["resegment", f"--id={RUN_ID}", f"--xenium-bundle={bundle}"]
        inputs = sorted(p for p in bundle.rglob("*") if p.is_file())
        if "expansion_um" in manifest:
            expansion = manifest["expansion_um"]
            if type(expansion) is not int or expansion < 0:
                raise ValueError("expansion_um must be a nonnegative integer")
            command.append(f"--expansion-distance={expansion}")
        reference_id = None
    else:
        engine = "spaceranger"
        prefix = required_text(manifest, "fastq_sample")
        fastqs = asset(root, manifest, "fastqs", directory=True)
        reference = asset(root, manifest, "reference", directory=True)
        reference_id = required_text(manifest, "reference_id")
        if not (reference / "reference.json").is_file():
            raise ValueError("Expected a 10x transcriptome reference with reference.json")
        inputs = fastq_pairs(fastqs, prefix) + sorted(p for p in reference.rglob("*") if p.is_file())
        command = ["count", f"--id={RUN_ID}", f"--fastqs={fastqs}", f"--sample={prefix}",
                   f"--transcriptome={reference}", "--create-bam=false", POLICIES[assay][policy]]
        for key, flag, required in (("image", "image", True), ("cytaimage", "cytaimage", assay == "probe" or platform == "Visium HD"),
                                    ("probe_set", "probe-set", assay == "probe"), ("slidefile", "slidefile", False)):
            if required or key in manifest:
                if key == "probe_set" and assay != "probe":
                    raise ValueError("Use probe_set only for probe-based assays")
                path = asset(root, manifest, key)
                command.append(f"--{flag}={path}")
                inputs.append(path)
        slide, area = required_text(manifest, "slide"), required_text(manifest, "area")
        command.extend((f"--slide={slide}", f"--area={area}"))
    for value, flag in ((cores, "localcores"), (memory_gb, "localmem")):
        if value is not None:
            if value < 1:
                raise ValueError(f"{flag} must be positive")
            command.append(f"--{flag}={value}")
    return sample, {"engine":engine, "expected_version":version, "assay":assay, "count_policy":policy,
                    "reference_id":reference_id, "command":[executable or engine, *command],
                    "source_record":link_source(manifest, root, sample)}, inputs


def output_assets(outs, platform):
    if platform == "Xenium":
        yield None, {"h5":outs / "cell_feature_matrix.h5", "cells":outs / "cells.csv.gz",
                     "experiment":outs / "experiment.xenium", "image":focus_image(outs)}
        return
    for size in ((2, 8, 16) if platform == "Visium HD" else (None,)):
        base = outs if size is None else outs / "binned_outputs" / f"square_{size:03}um"
        spatial = base / "spatial"
        names = ["tissue_positions.parquet"] if size else ["tissue_positions.csv", "tissue_positions_list.csv"]
        positions = next((spatial / name for name in names if (spatial / name).is_file()), spatial / names[0])
        yield size, {"h5":base / "filtered_feature_bc_matrix.h5", "positions":positions,
                     "scalefactors":spatial / "scalefactors_json.json", "image":spatial / "tissue_lowres_image.png"}


def save_json(path, value):
    pending = path.with_suffix(".tmp")
    pending.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    pending.replace(path)


def read_plan(manifest_path, executable, cores, memory_gb):
    manifest_path = manifest_path.resolve(strict=True)
    body = manifest_path.read_bytes()
    manifest = json.loads(body)
    if not isinstance(manifest, dict):
        raise ValueError("Expected a JSON object describing one sample and count policy")
    sample, plan, inputs = create_plan(manifest, manifest_path.parent, executable=executable,
                                       cores=cores, memory_gb=memory_gb)
    return body, sample, plan, inputs


def recount(manifest_path, out, *, executable=None, cores=None, memory_gb=None, plan_only=False):
    if plan_only:
        _, sample, plan, inputs = read_plan(manifest_path, executable, cores, memory_gb)
        return {"sample":sample, **plan, "input_files":[str(p) for p in inputs]}
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    journal = out / "recount.json"
    record = {"manifest":str(manifest_path), "stage":"preflight", "status":"running"}
    save_json(journal, record)
    try:
        body, sample, plan, inputs = read_plan(manifest_path, executable, cores, memory_gb)
        record.update(sample=sample, **plan)
        (out / "input_manifest.json").write_bytes(body)
        command = plan["command"]
        engine = shutil.which(command[0])
        if engine is None:
            raise ValueError(f"Required executable unavailable: {command[0]}")
        engine = command[0] = str(Path(engine).resolve(strict=True))
        version = subprocess.run([engine, "--version"], check=True, capture_output=True, text=True).stdout.strip()
        if not re.search(rf"(?<![\d.]){re.escape(plan['expected_version'])}(?![\w.+-])", version):
            raise ValueError(f"Engine version differs from manifest: {version!r}")
        plan["engine_version"] = version
        # Large input digests are stored once, outside the per-bin AnnData metadata.
        manifest = out / "inputs.json"
        save_json(manifest, {"input_manifest_sha256":sha256(out / "input_manifest.json"),
                            "files":[{"path":str(p), "bytes":p.stat().st_size, "sha256":sha256(p)} for p in inputs]})
        plan["inputs_sha256"] = sha256(manifest)
        record.update(plan, stage="counting")
        save_json(journal, record)
        with (out / "engine.log").open("w", encoding="utf8") as log:
            subprocess.run(command, cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True)
        record.update(stage="validation")
        save_json(journal, record)
        rows = []
        for size, paths in output_assets(out / RUN_ID / "outs", sample["platform"]):
            row = {**sample, "source":"local_recount", "status":"recount_ready_for_validation",
                   "recount":plan, "files":{role:{"path":str(p.relative_to(out)), "sha256":sha256(p)}
                                            for role, p in paths.items() if p.is_file()}}
            if size is not None:
                row.update(bin_um=size, sample_id=f"{sample['sample_id']}--{size}um")
            rows.append(row)
        report = out / "spatial_assets.json"
        save_json(report, {"request":{"organism":sample["organism"]}, "samples":rows})
        result = build([report], out / "dataset", local_recount=True)
        record["totals"] = result["totals"]
        if record["totals"]["built"] != len(rows):
            raise ValueError("One or more spatial outputs failed validation. See dataset/catalog.json")
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError, KeyboardInterrupt) as exc:
        record.update(status="failed", error=str(exc))
        save_json(journal, record)
        raise
    record.update(stage="complete", status="completed")
    save_json(journal, record)
    return record
