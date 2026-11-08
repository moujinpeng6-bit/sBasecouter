"""Prepare raw recount inputs from exact source matches and pinned public files."""

import hashlib
from http.client import HTTPException
import json
import os
from pathlib import Path, PurePosixPath
import re
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

from .build import HUMAN, sha256
from .checkpoint import output_lock
from .recount import read_plan, save_json


PUBLIC_HOSTS = {"cf.10xgenomics.com", "cdn.10xgenomics.com", "ftp.ncbi.nlm.nih.gov",
                "ftp.sra.ebi.ac.uk", "ftp.ebi.ac.uk", "ftp-trace.ncbi.nlm.nih.gov"}


def public_url(url):
    parts = urlsplit(url)
    mirror = parts.netloc == "s3-us-west-2.amazonaws.com" and parts.path.startswith("/10x.files/samples/")
    if parts.scheme != "https" or (parts.netloc not in PUBLIC_HOSTS and not mirror) or parts.fragment:
        raise ValueError("Expected a direct HTTPS URL on a supported public data host")
    return url


class PublicRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return super().redirect_request(req, fp, code, msg, headers, public_url(newurl))


def fetch_asset(spec, path):
    url = public_url(spec["url"])
    size, digest = spec["bytes"], spec["sha256"]
    if type(size) is not int or size <= 0 or not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError("A raw download must pin its byte length and SHA-256")
    if path.exists():
        if path.stat().st_size != size or sha256(path) != digest:
            raise ValueError(f"Existing file differs from the published recipe: {path.name}")
        return "reused"
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".partial")
    if part.is_symlink():
        raise ValueError("Partial downloads must be regular files")
    start = part.stat().st_size if part.is_file() else 0
    if start > size or (start == size and sha256(part) != digest):
        start = 0
    if start < size:
        headers = {"Accept-Encoding": "identity", "User-Agent": "scRecouter-spatial public data client"}
        if start:
            headers["Range"] = f"bytes={start}-"
        with build_opener(PublicRedirect()).open(Request(url, headers=headers), timeout=45) as response:
            if response.status == 206:
                if response.headers.get("Content-Range") != f"bytes {start}-{size - 1}/{size}":
                    raise OSError("Unexpected response range for the pinned file")
            elif response.status == 200:
                start = 0
            else:
                raise OSError(f"Unexpected download response: {response.status}")
            with part.open("ab" if start else "wb") as stream:
                remaining = size - start
                while remaining:
                    chunk = response.read(min(remaining, 1024**2))
                    if not chunk:
                        raise OSError("Download ended before the pinned length")
                    stream.write(chunk)
                    remaining -= len(chunk)
                if response.read(1):
                    raise OSError("Download exceeded the pinned length")
                stream.flush()
                os.fsync(stream.fileno())
    if sha256(part) != digest:
        raise ValueError(f"Downloaded file failed SHA-256: {path.name}")
    part.replace(path)
    return "downloaded"


def staged_path(root, name):
    relative = PurePosixPath(name)
    if not relative.parts or relative.is_absolute() or ".." in relative.parts or "\\" in name or ":" in name:
        raise ValueError("Download paths must stay within the sample's assets directory")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Download path resolves outside the sample's assets directory")
    return path


def prepare(report_path, recipe_path, out):
    report_path, recipe_path, out = report_path.resolve(), recipe_path.resolve(), out.resolve()
    report_bytes, recipe_bytes = report_path.read_bytes(), recipe_path.read_bytes()
    report = json.loads(report_bytes)
    recipes = json.loads(recipe_bytes)["samples"]
    index = {row["sample_id"]:row for row in recipes}
    if len(index) != len(recipes):
        raise ValueError("Each sample must have exactly one raw input recipe")
    rows = [row for row in report["samples"] if row.get("selected") is True or row.get("status") == "correspondence_passed"]
    if len({row["sample_id"] for row in rows}) != len(rows):
        raise ValueError("Selected report sample IDs must be unique")
    out.mkdir(parents=True, exist_ok=True)
    results = []
    with output_lock(out):
        identity = {"source_sha256":hashlib.sha256(report_bytes).hexdigest(),
                    "recipes_sha256":hashlib.sha256(recipe_bytes).hexdigest()}
        journal = out / "preparation.json"
        if journal.exists() and json.loads(journal.read_bytes()).get("inputs") != identity:
            raise ValueError("Choose a new preparation directory for changed source records or recipes")
        save_json(journal, {"inputs":identity, "samples":results})
        for row in rows:
            sample = row["sample_id"]
            record = {"sample_id":sample, "status":"raw_recipe_required"}
            recipe = index.get(sample)
            if recipe is not None:
                try:
                    record.update(prepare_sample(row, report, report_path, recipe, recipe_path, out))
                except (KeyError, TypeError, ValueError, OSError, HTTPException) as exc:
                    record.update(status="raw_preparation_failed", reason=str(exc))
            results.append(record)
            save_json(journal, {"inputs":identity, "samples":results})
    return results


def prepare_sample(row, report, report_path, recipe, recipe_path, out):
    keys = ("sample_id", "source_sample_id", "study_id", "platform")
    if any(not row.get(key) or row[key] != recipe.get(key) for key in keys):
        raise ValueError("Raw recipe does not match this study, sample and platform")
    declared = [value for value in (row.get("organism"), report.get("request", {}).get("organism"))
                if value and value != "*"]
    if not declared or any(value.casefold().strip() not in HUMAN for value in declared):
        raise ValueError("A human source must be confirmed before downloading raw assets")
    identity = hashlib.sha256(row["sample_id"].encode()).hexdigest()
    assets = (out / "assets" / identity).resolve()
    if not assets.is_relative_to(out):
        raise ValueError("Asset directory must stay under the preparation output")
    manifest_dir = (out / "manifests").resolve()
    if not manifest_dir.is_relative_to(out):
        raise ValueError("Manifest directory must stay under the preparation output")
    manifest_dir.mkdir(exist_ok=True)
    path = manifest_dir / f"{identity}.json"
    path.unlink(missing_ok=True)
    manifest = {**recipe["manifest"], **{key:row[key] for key in keys},
                "organism":declared[0], "source_url":row["source_url"],
                "source_report":Path(os.path.relpath(report_path, manifest_dir)).as_posix()}
    for key in ("fastqs", "image", "cytaimage", "probe_set", "slidefile", "loupe_alignment", "xenium_bundle"):
        if key in manifest:
            manifest[key] = Path(os.path.relpath(staged_path(assets, manifest[key]), manifest_dir)).as_posix()
    if "reference" in manifest:
        manifest["reference"] = Path(os.path.relpath((recipe_path.parent / manifest["reference"]).resolve(), manifest_dir)).as_posix()
    fetched = {name:fetch_asset(spec, staged_path(assets, name)) for name, spec in recipe["files"].items()}
    candidate = path.with_suffix(".pending")
    try:
        save_json(candidate, manifest)
        read_plan(candidate, None, None, None)
        candidate.replace(path)
    finally:
        candidate.unlink(missing_ok=True)
    return {"status":"ready", "manifest":str(path.relative_to(out)), "files":fetched}
