"""Extract only the public Visium files needed for the layer benchmark."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import file_digest
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
from zipfile import ZipFile

import pandas as pd

RECORD = "https://zenodo.org/records/6921620"
ARCHIVE = f"{RECORD}/files/snRNA-seq_and_spatial_transcriptomics.zip"
LABEL_SOURCE = "https://github.com/khodosevichlab/Schizophrenia20/blob/main/Notebook6_Visium.md#L537-L538"
LABEL_MEMBER = "visium_spatial_transcriptomics_histological_cortical_layer_annotations.zip"
COHORTS = {"control": ("MB7", "MB9", "MB11", "MB15", "MB19", "MB21", "MB18-2"),
           "schizophrenia": ("MB6", "MB12", "MB14", "MB22", "MB23", "MB50", "MB52")}
MEMBERS = {"h5": "filtered_feature_bc_matrix.h5", "positions": "spatial/tissue_positions_list.csv",
           "image": "spatial/tissue_lowres_image.png", "scalefactors": "spatial/scalefactors_json.json"}


def get_archive(output):
    from sragent_spatial.tenx_cache import RangeZip

    class Remote(RangeZip):
        def read(self, length=-1):
            if length < 0:
                length = self.size - self.tell()
            block = 4 * 1024**2
            return b"".join(super(Remote, self).read(min(block, length - i)) for i in range(0, length, block))

    return ZipFile(Remote(ARCHIVE, 14740207830, block_size=4 * 1024**2, max_bytes=96 * 1024**2,
                          cache=output / "ranges", timeout=45))


def prepare(output, selected, workers):
    from sragent_spatial.tenx_cache import verify_member

    def extract(source, member, local):
        info = source.getinfo(member)
        if verify_member(local, info):
            return
        with tempfile.NamedTemporaryFile(dir=local.parent, delete=False) as target:
            temp = Path(target.name)
            try:
                with source.open(info) as stream:
                    shutil.copyfileobj(stream, target)
            except Exception:
                target.close()
                temp.unlink(missing_ok=True)
                raise
        temp.replace(local)

    annotation_archive = output / "annotations.zip"
    with get_archive(output) as archive:
        extract(archive, LABEL_MEMBER, annotation_archive)

    def prepare_slice(entry):
        cohort, sample = entry
        destination = output / sample
        destination.mkdir(exist_ok=True)
        files, checksums = {}, {}
        with get_archive(output) as source:
            for role, member in MEMBERS.items():
                local = destination / Path(member).name
                extract(source, f"visium_spaceranger_output/{sample}/outs/{member}", local)
                files[role] = local.relative_to(output).as_posix()
                with local.open("rb") as stream:
                    checksums[role] = file_digest(stream, "sha256").hexdigest()
        label_file = destination / "labels.csv"
        with ZipFile(annotation_archive) as annotations:
            label_members = [name for name in annotations.namelist()
                             if re.fullmatch(re.escape(sample) + r"(?:_MB_\d{2}\.\d{2}\.\d{4})?\.csv", name)]
            if len(label_members) != 1:
                raise ValueError(f"Expected one author annotation table for {sample}")
            with annotations.open(label_members[0]) as stream:
                table = pd.read_csv(stream, dtype=str)
            if list(table.columns[:1]) != ["Barcode"] or len(table.columns) != 2:
                raise ValueError(f"Unexpected annotation format for {sample}")
            table.columns = ["barcode", "layer"]
            if not table["layer"].isin(["L1", "L2", "L3", "L4", "L5", "L6", "WM"]).all():
                raise ValueError(f"Unexpected annotation label for {sample}")
            table.to_csv(label_file, index=False)
        with label_file.open("rb") as stream:
            checksums["labels"] = file_digest(stream, "sha256").hexdigest()
        print(f"Prepared {sample}", flush=True)
        return {"slice": sample, "donor": sample, "cohort": cohort, "files": files,
                "labels": label_file.relative_to(output).as_posix(), "sha256": checksums}

    entries = [(cohort, sample) for cohort, samples in COHORTS.items() for sample in samples
               if not selected or sample in selected]
    if selected and selected != {sample for _, sample in entries}:
        raise ValueError("Unknown slice ID")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        samples = list(pool.map(prepare_slice, entries))
    return {"source": RECORD, "label_source": LABEL_SOURCE, "samples": samples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sragent-source", type=Path, required=True, help="Path to the SRAgent-spatial repository")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--sample", action="append", default=[])
    args = parser.parse_args()
    sys.path.insert(0, str((args.sragent_source / "src").resolve(strict=True)))
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    (root / "input.json").write_text(json.dumps(prepare(root, set(args.sample), args.workers), indent=2) + "\n", encoding="utf8")
