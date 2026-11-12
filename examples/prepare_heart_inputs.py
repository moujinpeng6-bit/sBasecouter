"""Resume the pinned 10x heart image and GRCh38 reference downloads."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import tarfile
import time

from tifffile import TiffFile

from discover_raw_visium import CDN, MIRROR, byte_range
from sbasecouter.checkpoint import output_lock
from sbasecouter.recount import save_json


PART = 8 * 1024**2
OBJECTS = {
    "he": ("samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Human_Heart/Visium_HD_6p5mm_Human_Heart_tissue_image.btf",
           5813436731, '"38795b628bce822b633467104d6f0f94-694"'),
    "reference": ("supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz",
                  11448662640, '"756a4eb76dccb0a7f9d909ace754c6ee-1365"'),
}


def download(out, role, workers):
    relative, size, etag = OBJECTS[role]
    identity = {"url": CDN + relative, "mirror": MIRROR + relative, "bytes": size, "etag": etag}
    target = out / Path(relative).name
    cache = out / "heart_segments" / role
    cache.mkdir(parents=True, exist_ok=True)
    starts = range(0, size, PART)

    def segment(start):
        length = min(PART, size - start)
        path = cache / f"{start}-{start + length - 1}.range"
        if path.is_file():
            with path.open("rb") as handle:
                digest = handle.readline(128).strip()
                body = handle.read(length + 1)
            if len(body) == length and hashlib.sha256(body).hexdigest().encode() == digest:
                return path
        for attempt in range(4):
            try:
                _, body = byte_range(identity["url"], start, length, identity)
                break
            except OSError:
                if attempt == 3:
                    raise
                time.sleep(attempt + 1)
        temp = path.with_suffix(".tmp")
        with temp.open("wb") as handle:
            handle.write(hashlib.sha256(body).hexdigest().encode() + b"\n")
            handle.write(body)
        temp.replace(path)
        return path

    with output_lock(cache):
        if target.exists():
            source = target
        else:
            source = target.with_name(target.name + ".assembling")
            with ThreadPoolExecutor(max_workers=workers) as pool:
                parts = list(pool.map(segment, starts))
            with source.open("wb") as out_file:
                for path in parts:
                    with path.open("rb") as handle:
                        digest = handle.readline(128).strip()
                        body = handle.read(PART + 1)
                    if hashlib.sha256(body).hexdigest().encode() != digest:
                        raise ValueError(f"Saved range changed: {path.name}")
                    out_file.write(body)
        if source.stat().st_size != size:
            raise ValueError("Assembled file does not match the published size")
        sha, multipart = hashlib.sha256(), hashlib.md5()
        with source.open("rb") as handle:
            while body := handle.read(PART):
                sha.update(body)
                multipart.update(hashlib.md5(body).digest())
        if f'"{multipart.hexdigest()}-{len(starts)}"' != etag:
            raise ValueError("Assembled file does not match the published multipart ETag")
        if source != target:
            source.replace(target)
        if role == "he":
            with TiffFile(target) as image:
                contents = {"image_shape": list(image.series[0].shape), "is_bigtiff": image.is_bigtiff}
        else:
            required = {"reference.json", "fasta/genome.fa", "genes/genes.gtf"}
            with tarfile.open(target, mode="r|gz") as archive:
                for entry in archive:
                    relative_name = entry.name.removeprefix("./").split("/", 1)[-1]
                    if entry.isfile() and entry.size > 0:
                        required.discard(relative_name)
                    if not required:
                        break
            if required:
                raise ValueError(f"Reference archive lacks required files: {sorted(required)}")
            contents = {"reference_id": "GRCh38-2020-A", "reference_files_verified": True}
        result = {"status": "complete", "asset": target.name, **identity, "sha256": sha.hexdigest(), **contents}
        save_json(out / f"{role}_windows_verified.json", result)
        return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--role", required=True, choices=OBJECTS)
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args()
    print(json.dumps(download(args.out.resolve(strict=True), args.role, args.workers)))
