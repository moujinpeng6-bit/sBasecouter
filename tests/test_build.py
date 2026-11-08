import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import anndata
import numpy as np
import pandas as pd
from PIL import Image

from sbasecouter.build import build, sha256


class BuildTests(unittest.TestCase):
    def fixture(self, root, platform="Visium", *, bin_um=None):
        barcodes = ["a", "b"] if platform == "Visium" else ["s_008um_00001_00001-1", "s_008um_00001_00002-1"]
        image = root/"tissue_hires_image.png"
        Image.new("L",(100,100)).save(image)
        assets = {"image":image, "positions":root/"positions.csv", "scalefactors":root/"scales.json",
                  "barcodes":root/"barcodes.tsv.gz", "features":root/"features.tsv.gz", "matrix":root/"matrix.mtx.gz"}
        text = {"barcodes":"\n".join(barcodes)+"\n",
                "features":"gene1\tGENE\tGene Expression\ncontrol1\tCONTROL\tNegative Control Probe\n",
                "matrix":"%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 4\n1 2 1\n2 1 2\n"}
        for role, content in text.items():
            with gzip.open(assets[role],"wt") as stream:
                stream.write(content)
        positions = [[barcodes[1],1,1,0,30,40], [barcodes[0],1,0,0,10,20]]
        if platform == "Visium":
            assets["positions"].write_text("".join(",".join(map(str,line))+"\n" for line in positions),encoding="utf8")
        else:
            assets["positions"] = root/"positions.parquet"
            columns = ["barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"]
            pd.DataFrame(positions,columns=columns).to_parquet(assets["positions"],index=False)
        assets["scalefactors"].write_text(json.dumps({"tissue_hires_scalef":2,"tissue_lowres_scalef":0.5}),encoding="utf8")
        row = {"sample_id":"a","platform":platform,"status":"correspondence_passed","bin_um":bin_um,
               "files":{role:{"path":path.name,"sha256":sha256(path)} for role,path in assets.items()}}
        report = root/"batch.json"
        report.write_text(json.dumps({"request":{"organism":"Homo sapiens"},"samples":[row]}),encoding="utf8")
        return report, row

    def test_counts_coordinate_order_and_hires_scale_are_preserved(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            report, _ = self.fixture(root)
            output = root/"out"
            catalog = build([report],output,normalize_total=100)
            result = catalog["samples"][0]
            self.assertEqual(result["status"],"built",result.get("reason"))
            data = anndata.read_h5ad(output/result["file"])
            self.assertEqual(data.obs_names.tolist(),["a","b"])
            np.testing.assert_array_equal(data.X.toarray(),[[4,2],[1,0]])
            np.testing.assert_array_equal(data.obsm["spatial"],[[20,10],[40,30]])
            np.testing.assert_array_equal(data.obsm["image_xy"],[[40,20],[80,60]])
            np.testing.assert_allclose(data.layers["log1p_total"].toarray(),[[np.log1p(100),0],[np.log1p(100),0]])

    def test_partial_hd_matrix_and_nonhuman_sample_are_skipped(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            report = root/"batch.json"
            report.write_text(json.dumps({"request":{"organism":"Homo sapiens"},"samples":[
                {"sample_id":"partial","platform":"Visium HD","status":"correspondence_passed","matrix_scope":"header_only"},
                {"sample_id":"mouse","platform":"Visium","status":"correspondence_passed","organism":"Mus musculus"},
            ]}),encoding="utf8")
            result = build([report],root/"out")
            self.assertEqual([r["reason"] for r in result["samples"]],
                             ["full_count_matrix_required","human_species_not_confirmed"])

    def test_rejects_changed_assets(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            report, _ = self.fixture(root)
            (root/"positions.csv").write_text("changed",encoding="utf8")
            result = build([report],root/"out")
            self.assertEqual(result["totals"]["failed"],1)
            self.assertIn("SHA-256",result["samples"][0]["reason"])

    def test_resume_adds_new_records_without_rewriting_verified_matrix(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            report, row = self.fixture(root)
            out = root / "out"
            first = build([report], out)["samples"][0]
            matrix = out / first["file"]
            before = matrix.stat().st_mtime_ns
            report.write_text(json.dumps({"request":{"organism":"Homo sapiens"},
                                          "samples":[row, {**row, "sample_id":"b"}]}),encoding="utf8")
            result = build([report], out, resume=True)
            self.assertEqual(result["totals"]["built"], 2)
            self.assertEqual(matrix.stat().st_mtime_ns, before)
            matrix.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "differs from the saved catalog"):
                build([report], out, resume=True)

    def test_complete_hd_bin_preserves_resolution_and_refuses_mismatch(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            report, row = self.fixture(root,platform="Visium HD",bin_um=8)
            sample = build([report],root/"valid")["samples"][0]
            self.assertEqual(sample["status"],"built",sample.get("reason"))
            self.assertEqual(sample["spatial_reference"]["bin_um"],8)
            row["bin_um"] = 16
            report.write_text(json.dumps({"request":{"organism":"Homo sapiens"},"samples":[row]}),encoding="utf8")
            sample = build([report],root/"invalid")["samples"][0]
            self.assertEqual(sample["status"],"failed")
            self.assertIn("bin size conflicts",sample["reason"])
