import json
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import anndata
import h5py
import numpy as np
import pandas as pd
from PIL import Image
from tifffile import imwrite

from sbasecouter.formats import VISIUM_FIELDS
from sbasecouter.recount import output_assets, recount


class RecountTests(unittest.TestCase):
    def manifest(self, root, platform, assay, policy):
        record = {"sample_id": "slice-a", "source_sample_id": "GSM123", "study_id": "GSE456",
                  "source_url": "https://example.org/slice-a", "organism": "Homo sapiens",
                  "platform": platform, "assay": assay, "count_policy": policy, "engine_version": "4.1.0"}
        if platform == "Xenium":
            bundle = root / "bundle"
            bundle.mkdir()
            for name in ("experiment.xenium", "transcripts.zarr.zip", "cells.zarr.zip", "morphology.ome.tif",
                         "morphology_focus.ome.tif", "metrics_summary.csv", "analysis_summary.html"):
                (bundle / name).write_text("fixture", encoding="utf8")
            record.update(xenium_bundle="bundle", expansion_um=3, engine_version="4.0.0")
        else:
            fastqs, reference = root / "reads", root / "reference"
            fastqs.mkdir()
            reference.mkdir()
            for lane in (1, 2):
                for read in (1, 2):
                    (fastqs / f"tissue_S1_L{lane:03}_R{read}_001.fastq.gz").write_bytes(b"fixture")
            (reference / "reference.json").write_text("{}", encoding="utf8")
            (root / "image.tif").write_bytes(b"fixture")
            record.update(fastqs="reads", fastq_sample="tissue", reference="reference", reference_id="GRCh38-demo",
                          image="image.tif", slide="V11-demo", area="A1")
            if assay == "probe" or platform == "Visium HD":
                (root / "cyta.tif").write_bytes(b"fixture")
                record["cytaimage"] = "cyta.tif"
            if assay == "probe":
                (root / "probes.csv").write_bytes(b"fixture")
                record["probe_set"] = "probes.csv"
        path = root / "sample.json"
        path.write_text(json.dumps(record), encoding="utf8")
        return path, record

    def write_engine_output(self, command, *, platform, version, omit=None, **kwargs):
        if command[1] == "--version":
            return subprocess.CompletedProcess(command, 0, stdout=f"vendor-{version}\n")
        outs = kwargs["cwd"] / "spatial_count" / "outs"
        for size, paths in output_assets(outs, platform):
            if size == omit and omit is not None:
                continue
            barcode = f"s_{size:03}um_00001_00001-1" if size else "cell-a"
            for path in paths.values():
                path.parent.mkdir(parents=True, exist_ok=True)
            with h5py.File(paths["h5"], "w") as file:
                matrix = file.create_group("matrix")
                matrix["data"], matrix["indices"], matrix["indptr"], matrix["shape"] = [4], [0], [0, 1], [1, 1]
                matrix["barcodes"] = [barcode.encode()]
                features = matrix.create_group("features")
                features["id"], features["name"], features["feature_type"] = [b"gene1"], [b"GENE"], [b"Gene Expression"]
            if platform == "Xenium":
                pd.DataFrame([[barcode, 2, 3]], columns=["cell_id", "x_centroid", "y_centroid"]).to_csv(paths["cells"], index=False)
                paths["experiment"].write_text('{"pixel_size":0.5}', encoding="utf8")
                image = outs / "morphology_focus/ch0000_dapi.ome.tif"
                image.parent.mkdir()
                imwrite(image, np.zeros((20, 20), dtype=np.uint8))
            else:
                positions = pd.DataFrame([[barcode, 1, 0, 0, 2, 3]], columns=VISIUM_FIELDS)
                if size:
                    positions.to_parquet(paths["positions"], index=False)
                else:
                    positions.to_csv(paths["positions"], index=False)
                paths["scalefactors"].write_text('{"tissue_lowres_scalef":0.5}', encoding="utf8")
                Image.new("L", (20, 20)).save(paths["image"])
        return subprocess.CompletedProcess(command, 0)

    def run_fixture(self, manifest, destination, *, platform, version, omit=None):
        with patch("sbasecouter.recount.shutil.which", return_value=str(manifest)), patch(
            "sbasecouter.recount.subprocess.run",
            side_effect=lambda command, **kw: self.write_engine_output(command, platform=platform, version=version, omit=omit, **kw),
        ):
            return recount(manifest, destination)

    def test_all_reads_and_count_policy_with_exact_source_record(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, manifest = self.manifest(root, "Visium", "three_prime", "exonic")
            row = {**manifest, "observed_sample_state": "healthy"}
            (root / "batch.json").write_text(json.dumps({"samples": [row]}), encoding="utf8")
            manifest["source_report"] = "batch.json"
            path.write_text(json.dumps(manifest), encoding="utf8")
            plan = recount(path, root / "plan-only", plan_only=True)
            self.assertIn("--include-introns=false", plan["command"])
            self.assertEqual(sum(p.endswith(".fastq.gz") for p in plan["input_files"]), 4)
            self.assertEqual(plan["sample"]["observed_sample_state"], "healthy")
            self.assertEqual(plan["source_record"]["link"], "sragent_spatial_report")
            result = self.run_fixture(path, root / "valid", platform="Visium", version="4.1.0")
            self.assertEqual(result["totals"]["built"], 1)
            dataset = root / "valid/dataset"
            data = anndata.read_h5ad(next(dataset.glob("*.h5ad")))
            self.assertEqual(data.X[0, 0], 4)
            self.assertEqual(json.loads(data.uns["source_record_json"])["sample"]["observed_sample_state"], "healthy")
            manifest["count_policy"] = "exonic_intronic"
            path.write_text(json.dumps(manifest), encoding="utf8")
            self.assertIn("--include-introns=true", recount(path, root / "plan-only", plan_only=True)["command"])
            (root / "reads/tissue_S1_L002_R1_001.fastq.gz").unlink()
            with self.assertRaisesRegex(ValueError, "R1/R2"):
                recount(path, root / "bad")
            self.assertEqual(json.loads((root / "bad/recount.json").read_text())["stage"], "preflight")

    def test_hd_outputs_and_missing_resolution_remain_distinct(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, _ = self.manifest(root, "Visium HD", "probe", "all_nondeprecated_probes")
            result = self.run_fixture(manifest, root / "valid", platform="Visium HD", version="4.1.0")
            self.assertEqual(result["totals"]["built"], 3)
            self.assertIn("--filter-probes=false", result["command"])
            catalog = json.loads((root / "valid/dataset/catalog.json").read_text())
            self.assertEqual([r["spatial_reference"]["bin_um"] for r in catalog["samples"]], [2, 8, 16])
            with self.assertRaisesRegex(ValueError, "failed validation"):
                self.run_fixture(manifest, root / "partial", platform="Visium HD", version="4.1.0", omit=16)
            report = json.loads((root / "partial/recount.json").read_text())
            self.assertEqual((report["status"], report["totals"]["built"], report["totals"]["failed"]), ("failed", 2, 1))
            resolution = Path("spatial_count/outs/binned_outputs/square_016um")
            shutil.copytree(root / "valid" / resolution, root / "partial" / resolution)
            with patch("sbasecouter.recount.shutil.which", return_value=str(manifest)), patch(
                "sbasecouter.recount.subprocess.run", return_value=subprocess.CompletedProcess([], 0, stdout="vendor-4.1.0"),
            ) as vendor:
                recovered = recount(manifest, root / "partial", resume=True)
            self.assertEqual(recovered["totals"]["built"], 3)
            self.assertEqual(vendor.call_count, 1)

    def test_xenium_output_has_real_count_provenance_and_spatial_coordinates(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, _ = self.manifest(root, "Xenium", "in_situ", "resegment")
            with self.assertRaisesRegex(ValueError, "version differs"):
                self.run_fixture(manifest, root / "wrong-version", platform="Xenium", version="4.0.1")
            result = self.run_fixture(manifest, root / "valid", platform="Xenium", version="4.0.0")
            self.assertIn("--expansion-distance=3", result["command"])
            dataset = root / "valid/dataset"
            data = anndata.read_h5ad(next(dataset.glob("*.h5ad")))
            self.assertEqual(data.obs_names.tolist(), ["cell-a"])
            np.testing.assert_array_equal(data.obsm["image_xy"], [[4, 6]])
            source = json.loads(data.uns["source_record_json"])
            self.assertEqual(source["matrix_kind"], "local_recount")
            self.assertEqual(source["recount"]["count_policy"], "resegment")

    def test_resume_reuses_engine_outputs_and_rejects_changed_fastq(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest, _ = self.manifest(root, "Visium", "three_prime", "exonic")
            out = root / "out"
            with patch("sbasecouter.recount.build", side_effect=KeyboardInterrupt):
                with self.assertRaises(KeyboardInterrupt):
                    self.run_fixture(manifest, out, platform="Visium", version="4.1.0")
            self.assertTrue(json.loads((out / "recount.json").read_text())["engine_completed"])
            with patch("sbasecouter.recount.shutil.which", return_value=str(manifest)), patch(
                "sbasecouter.recount.subprocess.run", return_value=subprocess.CompletedProcess([], 0, stdout="vendor-4.1.0"),
            ) as vendor:
                self.assertEqual(recount(manifest, out, resume=True)["status"], "completed")
                self.assertEqual(recount(manifest, out, resume=True)["totals"]["built"], 1)
                self.assertTrue(all(call.args[0][1] == "--version" for call in vendor.call_args_list))
                next((root / "reads").glob("*.fastq.gz")).write_bytes(b"changed")
                with self.assertRaisesRegex(ValueError, "different inputs"):
                    recount(manifest, out, resume=True)
            self.assertEqual(json.loads((out / "recount.json").read_text())["status"], "completed")
