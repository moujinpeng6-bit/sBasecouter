"""Resolve public 10x Visium HD raw inputs from SRAgent-spatial searches."""

import argparse
import html
import io
import json
from pathlib import Path
import re
import struct
import sys
import tarfile
from urllib.error import HTTPError
from urllib.request import Request, urlopen


CDN = "https://cf.10xgenomics.com/"
MIRROR = "https://s3.us-west-2.amazonaws.com/10x.files/"
SUMMARY = re.compile(r"samples/spatial-exp/(\d+\.\d+\.\d+)/([\w-]+)/\2_web_summary\.html")


def byte_range(url, start=0, length=1, expected=None):
    if not url.startswith(CDN):
        raise ValueError("Expected an official 10x file URL")
    span = f"{start}-{start + length - 1}"
    headers = {"Range": f"bytes={span}", "Accept-Encoding": "identity"}
    if expected:
        headers["If-Match"] = expected["etag"]
    source = MIRROR + url.removeprefix(CDN)
    with urlopen(Request(f"{source}?sragent_range={span}", headers=headers), timeout=25) as response:
        match = re.fullmatch(rf"bytes {span}/(\d+)", response.headers.get("Content-Range", ""))
        etag = response.headers.get("ETag", "")
        if response.status != 206 or not match or not re.fullmatch(r'"[a-f0-9]+(?:-\d+)?"', etag):
            raise ValueError("Range response did not pin the requested source object")
        identity = {"url": url, "mirror": source, "bytes": int(match[1]), "etag": etag}
        data = response.read(length)
    if len(data) != length or (expected and identity != expected):
        raise ValueError("Source size, range or ETag changed")
    return identity, data


def cytassist_tags(asset):
    url = asset["url"]
    _, header = byte_range(url, length=8, expected=asset)
    endian = {b"II": "<", b"MM": ">"}[header[:2]]
    marker, offset = struct.unpack(endian + "HI", header[2:])
    if marker != 42:
        raise ValueError("Expected the publisher's CytAssist TIFF")
    _, raw = byte_range(url, offset, 2, asset)
    count, = struct.unpack(endian + "H", raw)
    if not 1 <= count <= 512:
        raise ValueError("Unexpected CytAssist TIFF tag count")
    _, table = byte_range(url, offset + 2, count * 12, asset)
    tags = {}
    for entry in struct.iter_unpack(endian + "HHI4s", table):
        tag, kind, size, inline = entry
        if tag in (65002, 65012) and kind == 2 and 0 < size <= 4096:
            start, = struct.unpack(endian + "I", inline)
            content = inline[:size] if size <= 4 else byte_range(url, start, size, asset)[1]
            tags[tag] = content.decode("ascii").rstrip("\0")
    return tags


def resolve(card, result, slide, area):
    summary = result["url"]
    match = SUMMARY.fullmatch(summary.removeprefix(CDN)) if summary.startswith(CDN) else None
    if not match:
        raise ValueError("Expected an official, sample-specific Space Ranger summary")
    version, name = match.groups()
    command = next((line for line in html.unescape(result["content"]).splitlines() if "--fastqs=" in line), "")
    args = dict(re.findall(r"--([a-z][a-z0-9-]*)=([^\s\"'`]+)", command))
    if args.get("slide", "").casefold() != slide.casefold() or args.get("area", "").replace("-", "").upper() != area:
        raise ValueError("Publisher card and summary disagree on slide or capture area")
    if not all(args.get(flag) for flag in ("fastqs", "image", "cytaimage", "transcriptome")):
        raise ValueError("Publisher run summary lacks required input flags")
    reference = Path(args["transcriptome"]).name
    if not re.fullmatch(r"refdata-gex-GRCh38-\d{4}-[A-Z]", reference):
        raise ValueError("Expected a named human 10x reference")
    prefix = summary.rsplit("/", 1)[0] + "/"
    names = {"fastq_archive": [f"{name}_fastqs.tar"], "cytassist": [f"{name}_image.tif", Path(args["cytaimage"]).name],
             "microscope": [f"{name}_tissue_image.btf", Path(args["image"]).name]}
    if "probe-set" in args:
        names["probe_set"] = [f"{name}_probe_set.csv", Path(args["probe-set"]).name]
    assets = {}
    for role, filenames in names.items():
        for filename in dict.fromkeys(filenames):
            if not re.fullmatch(r"[A-Za-z0-9._-]+", filename) or not filename.startswith((name + "_", "Visium_Human_Transcriptome_Probe_Set_")):
                continue
            try:
                assets[role] = byte_range(prefix + filename)[0]
                break
            except HTTPError as exc:
                if exc.code not in {403, 404}:
                    raise
        if role not in assets:
            raise ValueError(f"No readable publisher object for {role}")
    assets["reference"] = byte_range(f"{CDN}supp/cell-exp/{reference}.tar.gz")[0]
    fastq = assets["fastq_archive"]
    _, head = byte_range(fastq["url"], length=16384, expected=fastq)
    sample = args.get("sample", name)
    with tarfile.open(fileobj=io.BytesIO(head), mode="r|") as archive:
        first = next((member.name for member in archive if member.isfile() and member.name.endswith(".fastq.gz")), "")
    if not first or not re.fullmatch(re.escape(sample) + r"_S\d+(?:_L\d{3})?_[RI][12]_\d{3}\.fastq\.gz", first.rsplit("/", 1)[-1]):
        raise ValueError("FASTQ archive member disagrees with the summary sample prefix")
    tags = cytassist_tags(assets["cytassist"])
    if slide.casefold() not in tags.get(65002, "").casefold() or "cytassist" not in tags.get(65012, "").casefold():
        raise ValueError("CytAssist TIFF tags disagree with the publisher card")
    return {"status": "public_sources_verified", "sample_id": card["sample_id"], "card": card["source_url"],
            "publisher_summary": summary, "organism": "Homo sapiens", "platform": "Visium HD", "assay": "probe" if "probe-set" in args else "three_prime",
            "slide": slide, "area": area, "engine_version": version, "reference_id": reference.removeprefix("refdata-gex-"),
            "fastq_sample": sample, "publisher_fastq_directories": [Path(p).name for p in args["fastqs"].split(",")],
            "first_fastq_member_checked": first, "cytassist_slide_tag_verified": True, "assets": assets}


def discover(tissue, search, card_row, keywords="", search_label=""):
    card_query = f"site:10xgenomics.com/datasets Human Visium HD {tissue} {keywords}".rstrip()
    summary_query = f'site:cf.10xgenomics.com "Visium HD" "Human {tissue}" web_summary {keywords}'.rstrip()
    cards = search(card_query, 12)
    summaries = search(summary_query, 12)
    queries = [card_query, summary_query]
    if search_label:
        public_query = f'site:10xgenomics.com "{search_label}" "{keywords}" "Visium HD"'
        queries.append(public_query)
        related = search(public_query, 12)
        cards = related + cards
        summaries = related + summaries
    outcomes = []
    seen_cards = set()
    for candidate in cards:
        row = card_row(candidate, tissue, "Visium HD")
        if row is None or row["sample_id"] in seen_cards:
            continue
        seen_cards.add(row["sample_id"])
        body = candidate.get("content", "")
        slides = set(re.findall(r"Slide serial number:\s*([\w-]+)", body, re.I))
        areas = {value.replace("-", "").upper() for value in re.findall(r"\bArea:\s*([A-D](?:-?[1-4])?)\b", body, re.I)}
        if len(slides) != 1 or len(areas) != 1:
            outcomes.append({"card": row["source_url"], "status": "sample_identity_unresolved"})
            continue
        slide, = slides
        area, = areas
        pattern = re.compile(r"--slide=" + re.escape(slide) + r"(?=\s|$)", re.I)
        matches = [summary for summary in summaries if pattern.search(summary.get("content", ""))]
        if not matches:
            targeted = f'site:cf.10xgenomics.com "{slide}" web_summary'
            queries.append(targeted)
            matches = [summary for summary in search(targeted, 6) if pattern.search(summary.get("content", ""))]
        if not matches:
            outcomes.append({"card": row["source_url"], "status": "publisher_summary_unresolved"})
        official = {item["url"]: item for item in matches if item.get("url", "").startswith(CDN)
                    and SUMMARY.fullmatch(item["url"].removeprefix(CDN))}
        explicit = {url: item for url, item in official.items() if url in body}
        if len(explicit or official) > 1:
            outcomes.append({"card": row["source_url"], "status": "ambiguous_publisher_summary",
                             "summaries": sorted(explicit or official)})
            continue
        for summary in (explicit or official).values():
            try:
                return {"queries": queries, "manifest": resolve(row, summary, slide, area), "other_candidates": outcomes}
            except (KeyError, OSError, ValueError, tarfile.TarError) as exc:
                outcomes.append({"card": row["source_url"], "status": "verification_failed",
                                 "reason": f"{type(exc).__name__}: {exc}"})
    return {"queries": queries, "manifest": None, "other_candidates": outcomes}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sragent-source", type=Path, required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--keywords", default="")
    parser.add_argument("--search-label", default="", help="Optional sample description used to refine the public index search")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str((args.sragent_source / "src").resolve(strict=True)))
    from sragent_spatial.tenx_index import card_row, lookup_public_cards
    report = discover(args.tissue, lookup_public_cards, card_row, args.keywords, args.search_label)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")
    sys.exit(0 if report["manifest"] else 1)
