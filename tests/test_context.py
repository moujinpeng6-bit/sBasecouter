from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

from sbasecouter.context import build_context, contextualize, neighbor_graph


class ContextTests(unittest.TestCase):
    def fixture(self):
        data = anndata.AnnData(
            sparse.csr_matrix([[2,0,4],[0,5,0],[3,0,0],[0,0,1],[0,3,0]]),
            obs=pd.DataFrame(index=["a","b","c","zero","isolated"]),
            var=pd.DataFrame({"feature_type":["Gene Expression","Gene Expression","Negative Control Probe"]},
                             index=["g1","g2","negative"]))
        data.obsm["spatial"] = np.array([[0,0],[1,0],[2,0],[3,0],[20,0]],dtype=float)
        data.uns["spatial_reference"] = {"observation_unit":"cell","native_unit":"micrometre","slice_id":"section-a"}
        return data

    def test_real_entry_point_preserves_counts_and_aligns_cell_types(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = self.fixture()
            labels = root/"labels.csv"
            labels.write_text("slice_id,barcode,cell_type\nsection-a,c,B\nsection-a,a,T\nsection-a,b,T\n",encoding="utf8")
            source, output = root/"input.h5ad", root/"output.h5ad"
            data.write_h5ad(source)
            result = build_context(source, output, neighbors=8, radius=3, alpha=0.25, same_type_weight=2,
                                   label_file=labels, label_source="example annotation")
            reopened = anndata.read_h5ad(output)
            self.assertEqual((reopened.X != data.X).nnz,0)
            self.assertTrue(reopened.obs_names.equals(data.obs_names))
            self.assertTrue(reopened.obs["spatial_context_id"].is_unique)
            self.assertEqual(result["cells_with_known_types"],3)
            self.assertEqual(result["neighbor_only_estimates"],1)
            self.assertEqual(reopened.uns["spatial_context"]["labels_reference"], "example annotation")
            self.assertEqual(reopened.uns["spatial_context"]["neighbor_type_columns"].tolist(),["B","T","unknown"])
            np.testing.assert_allclose(reopened.obsm["neighbor_type_fractions"].toarray()[0],[1/3,1/3,1/3])
            self.assertEqual(reopened.obsp["spatial_neighbors"].getrow(0).nnz,3)
            self.assertEqual(reopened.obsp["spatial_rna_weights"].getrow(0).nnz,2)
            layer = reopened.layers["spatial_context_expression"].toarray()
            np.testing.assert_allclose(layer[0],[25000/3,5000/3,0],rtol=1e-6)
            np.testing.assert_allclose(layer[4],[0,10000,0])

    def test_no_labels_stays_unknown_and_cross_slice_labels_fail(self):
        data = self.fixture()
        result = contextualize(data,neighbors=8,radius=3,alpha=0.25,same_type_weight=2)
        self.assertEqual(result["cells_with_known_types"],0)
        self.assertTrue(data.obs["cell_type"].eq("unknown").all())
        with TemporaryDirectory() as tmp:
            labels = Path(tmp)/"labels.csv"
            labels.write_text("slice_id,barcode,cell_type\nsection-b,a,T\n",encoding="utf8")
            with self.assertRaisesRegex(ValueError,"slice ID"):
                contextualize(self.fixture(),neighbors=8,radius=3,alpha=0.25,same_type_weight=2,label_file=labels)

    def test_duplicate_coordinates_do_not_include_self_or_exceed_k(self):
        graph = neighbor_graph(np.zeros((8,2)),2,1)
        self.assertFalse(graph.diagonal().any())
        np.testing.assert_array_equal(np.diff(graph.indptr),np.full(8,2))

    def test_zero_blend_keeps_zero_rna_cells_unfilled(self):
        data = self.fixture()
        result = contextualize(data, neighbors=8, radius=3, alpha=0, same_type_weight=2)
        self.assertEqual(result["neighbor_only_estimates"], 0)
        self.assertEqual(data.layers["spatial_context_expression"].getrow(3).nnz, 0)
        self.assertGreater(data.layers["neighbor_expression"].getrow(3).nnz, 0)

    def test_refuses_spots_for_cell_neighborhoods(self):
        data = self.fixture()
        data.uns["spatial_reference"]["observation_unit"] = "spot"
        with self.assertRaisesRegex(ValueError,"cell-resolved"):
            contextualize(data,neighbors=8,radius=3,alpha=0.25,same_type_weight=2)
