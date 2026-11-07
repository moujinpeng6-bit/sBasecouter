"""10x expression formats and platform-specific spatial coordinates."""

import json
import math
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse
from scipy.io import mmread
from tifffile import TiffFile


def expression(files):
    if "h5" in files:
        with h5py.File(files["h5"]) as handle:
            source = handle["matrix"]
            data, indices, pointers = (source[name][:] for name in ("data", "indices", "indptr"))
            shape = tuple(map(int, source["shape"][:]))
            barcodes = source["barcodes"].asstr()[:]
            metadata = source["features"]
            features = pd.DataFrame({
                "feature_id": metadata["id"].asstr()[:],
                "feature_name": metadata["name"].asstr()[:],
                "feature_type": metadata["feature_type"].asstr()[:],
            })
            if "genome" in metadata:
                features["genome"] = metadata["genome"].asstr()[:]
            # In 10x HDF5, features × observations are stored as CSC. The same arrays
            # represent the transposed observations × features matrix as CSR.
            counts = sparse.csr_matrix((data, indices, pointers), shape=shape[::-1])
    else:
        counts = mmread(files["matrix"]).T
        if not sparse.issparse(counts):
            raise ValueError("Expected a sparse Matrix Market count matrix")
        barcodes = pd.read_csv(files["barcodes"], sep="\t", header=None, dtype=str, keep_default_na=False)[0]
        features = pd.read_csv(files["features"], sep="\t", header=None, dtype=str, keep_default_na=False)
        features = features.iloc[:, :3]
        if features.shape[1] == 2:
            features[2] = "Gene Expression"
        features.columns = ["feature_id", "feature_name", "feature_type"]
    barcodes = pd.Index(barcodes, dtype=str, name="barcode")
    if counts.shape != (len(barcodes), len(features)) or not all(counts.shape):
        raise ValueError("Matrix shape does not match barcodes and features")
    if not barcodes.is_unique or (barcodes == "").any() or (features["feature_id"] == "").any():
        raise ValueError("Empty or duplicate barcode, or missing feature ID")
    features.index = pd.Index(features["feature_id"], dtype=str, name="feature_key")
    if not features.index.is_unique:
        features.index = pd.Index(features["feature_type"] + ":" + features["feature_id"], name="feature_key")
    if not features.index.is_unique:
        raise ValueError("Duplicate feature ID within a feature type")
    data = counts.data
    if (data.dtype.kind not in "iuf" or not np.isfinite(data).all()
            or (data < 0).any() or (data > np.iinfo(np.int32).max).any()
            or (data.dtype.kind == "f" and (data != np.floor(data)).any())):
        raise ValueError("Expression must contain nonnegative integer counts in the supported range")
    counts = counts.astype(np.int32 if counts.has_canonical_format else np.int64, copy=False)
    counts = sparse.csr_matrix(counts)
    counts.check_format(full_check=True)
    counts.sum_duplicates()
    counts.eliminate_zeros()
    return counts, barcodes, features


VISIUM_FIELDS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]


def spatial(files, platform, barcodes, bin_um=None):
    image = files["image"]
    if platform == "Xenium":
        table = pd.read_csv(files["cells"], dtype={"cell_id":str}, keep_default_na=False)
        key, axes = "cell_id", ["x_centroid", "y_centroid"]
        metadata = json.loads(files["experiment"].read_text(encoding="utf-8-sig"))
        pixel_size = float(metadata["pixel_size"])
        if not math.isfinite(pixel_size) or pixel_size <= 0:
            raise ValueError("Invalid Xenium micrometres per pixel")
        scale = 1 / pixel_size
        with TiffFile(image) as tif:
            series = tif.series[0]
            width, height = (series.shape[series.axes.index(axis)] for axis in "XY")
        reference = {"observation_unit":"cell", "native_unit":"micrometre", "micrometres_per_image_pixel":pixel_size}
    else:
        path = files["positions"]
        if path.suffix == ".parquet":
            table = pd.read_parquet(path, columns=VISIUM_FIELDS)
        else:
            table = pd.read_csv(path, dtype={"barcode":str}, keep_default_na=False)
            if "barcode" not in table:
                table = pd.read_csv(path, names=VISIUM_FIELDS, header=None, dtype={"barcode":str}, keep_default_na=False)
        key, axes = "barcode", ["pxl_col_in_fullres", "pxl_row_in_fullres"]
        metadata = json.loads(files["scalefactors"].read_text(encoding="utf-8-sig"))
        image_name = Path(image).name.removesuffix(".gz")
        scale_key = {"tissue_lowres_image.png":"tissue_lowres_scalef",
                     "tissue_hires_image.png":"tissue_hires_scalef"}.get(image_name)
        if scale_key is None:
            raise ValueError("Expected the publisher's named Visium low-resolution or high-resolution PNG")
        scale = float(metadata[scale_key])
        with Image.open(image) as png:
            width, height = png.size
        reference = {"observation_unit":"spot", "native_unit":"full_resolution_image_pixel", "scale_factor_field":scale_key}
        if platform == "Visium HD":
            if path.suffix != ".parquet":
                raise ValueError("Expected the native Visium HD positions Parquet")
            sizes = barcodes.to_series().str.extract(r"^s\d*_(\d{3})um_\d+_\d+(?:-\d+)?$", expand=False)
            if sizes.isna().any() or sizes.nunique() != 1:
                raise ValueError("Visium HD barcodes must declare one native bin size")
            size = int(sizes.iloc[0])
            if size <= 0 or any(x is not None and float(x) != size for x in (bin_um, metadata.get("bin_size_um"))):
                raise ValueError("Visium HD bin size conflicts with report or scale factors")
            reference.update(observation_unit="bin", bin_um=size)
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid scale from native coordinates to image pixels")
    if table[key].isna().any() or (table[key] == "").any() or table[key].duplicated().any():
        raise ValueError("Empty or duplicate coordinate barcode")
    table = table.set_index(key)
    missing = barcodes.difference(table.index)
    if len(missing) or (platform == "Xenium" and len(table) != len(barcodes)):
        raise ValueError(f"Expression and coordinates disagree ({len(missing)} expression IDs lack coordinates)")
    joined = table.loc[barcodes]
    native_xy = joined[axes].to_numpy(dtype=np.float64)
    image_xy = native_xy * scale
    if not np.isfinite(image_xy).all():
        raise ValueError("Invalid coordinates after image scaling")
    obs = pd.DataFrame(index=barcodes)
    if platform != "Xenium":
        tissue = joined["in_tissue"].to_numpy(dtype=float)
        if not np.isin(tissue, [0, 1]).all():
            raise ValueError("Visium in_tissue must be 0 or 1")
        obs["in_tissue"] = tissue.astype(bool)
        for axis in ("array_row", "array_col"):
            grid = joined[axis].to_numpy(dtype=float)
            if not np.isfinite(grid).all() or (grid < 0).any() or (grid != np.floor(grid)).any():
                raise ValueError("Invalid Visium array index")
            obs[axis] = grid.astype(np.int32)
    checked = image_xy[obs["in_tissue"].to_numpy()] if "in_tissue" in obs else image_xy
    if not len(checked) or (checked < 0).any() or (checked >= [width, height]).any():
        raise ValueError("No tissue coordinates or mapped coordinates outside the image")
    reference.update(axis_order="x,y", image_width=int(width), image_height=int(height),
                     image_pixels_per_native_unit=scale, coordinates_checked=len(checked),
                     unused_coordinate_rows=len(table)-len(barcodes))
    return obs, native_xy, image_xy, reference
