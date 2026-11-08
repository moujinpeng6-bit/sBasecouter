import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from sbasecouter.prepare import fetch_asset, prepare, public_url, staged_path


class PrepareTests(unittest.TestCase):
    def test_partial_download_is_resumed_and_verified_before_reuse(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "asset"
            path.with_suffix(".partial").write_bytes(b"abc")
            spec = {"url":"https://cf.10xgenomics.com/samples/fixture", "bytes":6,
                    "sha256":hashlib.sha256(b"abcdef").hexdigest()}
            response = io.BytesIO(b"def")
            response.status = 206
            response.headers = {"Content-Range":"bytes 3-5/6"}
            with patch("sbasecouter.prepare.build_opener") as opener:
                opener.return_value.open.return_value = response
                self.assertEqual(fetch_asset(spec,path), "downloaded")
                self.assertEqual(opener.return_value.open.call_args.args[0].get_header("Range"), "bytes=3-")
                self.assertEqual(fetch_asset(spec,path), "reused")
                self.assertEqual(opener.return_value.open.call_count, 1)
            path.write_bytes(b"broken")
            with self.assertRaisesRegex(ValueError,"differs"):
                fetch_asset(spec,path)

    def test_untrusted_url_and_path_are_rejected(self):
        with self.assertRaises(ValueError):
            public_url("https://cf.10xgenomics.com@example.org/fixture")
        with self.assertRaises(ValueError):
            staged_path(Path.cwd(), "../fixture")

    def test_missing_and_mismatched_recipe_are_reported_per_sample(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            report, recipes = root / "batch.json", root / "recipes.json"
            row = {"sample_id":"x", "study_id":"study", "source_sample_id":"id", "platform":"Visium", "selected":True}
            report.write_text(json.dumps({"request":{"organism":"Homo sapiens"},"samples":[row,{**row,"sample_id":"y"}]}))
            recipes.write_text(json.dumps({"samples":[{**row,"study_id":"different"}]}))
            with patch("sbasecouter.prepare.fetch_asset") as fetch:
                results = prepare(report,recipes,root / "prepared")
            self.assertEqual([row["status"] for row in results], ["raw_preparation_failed", "raw_recipe_required"])
            self.assertIn("does not match", results[0]["reason"])
            fetch.assert_not_called()

    def test_ready_manifest_resolves_raw_assets_and_source(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            report, recipes = root / "batch.json", root / "recipes.json"
            row = {"sample_id":"x", "study_id":"study", "source_sample_id":"id", "platform":"Visium",
                   "selected":True,"source_url":"https://www.10xgenomics.com/datasets/example"}
            report.write_text(json.dumps({"request":{"organism":"Homo sapiens"},"samples":[row]}))
            (root / "reference").mkdir()
            (root / "reference/reference.json").write_text("{}")
            assets = ["reads/test_S1_R1_001.fastq.gz", "reads/test_S1_R2_001.fastq.gz", "image.tif"]
            manifest = {"assay":"three_prime", "count_policy":"exonic", "engine_version":"4.1.0",
                        "fastqs":"reads", "fastq_sample":"test", "image":"image.tif", "reference":"reference",
                        "reference_id":"GRCh38-example", "slide":"V11-example", "area":"A1"}
            recipes.write_text(json.dumps({"samples":[{**row,"manifest":manifest,"files":dict.fromkeys(assets,{})}]}))

            def download(_spec, path):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"fixture")
                return "downloaded"

            out = root / "prepared"
            with patch("sbasecouter.prepare.fetch_asset", side_effect=download):
                record, = prepare(report, recipes, out)
            self.assertEqual(record["status"], "ready")
            self.assertTrue((out / record["manifest"]).is_file())
