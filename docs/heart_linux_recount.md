# 人类心脏 Visium HD 的 Linux 重计数准备

样本是 10x 公布的人类健康心脏 FFPE 切片。原始运行使用 Space Ranger 4.1.0、参考库 `GRCh38-2020-A`、人类 Visium 探针集 v2.0、载玻片 `H1-ZX2RP3Y` 和捕获区 `A1`。[样本资料](https://www.10xgenomics.com/cn/datasets/visium-hd-cytassist-6p5mm-human-heart) · [原始运行摘要](https://cf.10xgenomics.com/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Human_Heart/Visium_HD_6p5mm_Human_Heart_web_summary.html)

## 服务器和输入

- 准备 64 位 Linux，例如 Ubuntu 22.04 LTS，采用支持 AVX 的 Intel 或 AMD CPU。10x 的最低要求为 8 核、64 GB 内存和 1 TB 可用磁盘。对这张全量切片，按官方推荐准备 32 核、128 GB 内存，并给工作目录留至少 1 TB 可用空间。[10x 系统要求](https://www.10xgenomics.com/support/jp/software/space-ranger/downloads/space-ranger-system-requirements)
- 安装 Space Ranger **4.1.0** 和 Python 3.11 或更高版本。将 `refdata-gex-GRCh38-2020-A` 参考库解包成目录，目录内应有 `reference.json`。[Space Ranger 与参考库下载中心](https://www.10xgenomics.com/support/software/space-ranger/downloads)
- 准备 `Visium_HD_6p5mm_Human_Heart_fastqs.tar`。解包后以 `Visium_HD_6p5mm_Human_Heart_fastqs` 为 FASTQ 目录，其中有 L001、L002 两条 lane，每条 lane 有 R1、R2 及 I1、I2。样本名前缀为 `Visium_HD_6p5mm_Human_Heart`。[10x 的 FASTQ 输入约定](https://www.10xgenomics.com/support/software/space-ranger/latest/tutorials/inputs/fastqs-specifying-fastqs)
- 将 `Visium_HD_6p5mm_Human_Heart_image.tif` 用作 **CytAssist 图像**。文件与官方空间归档中的 `spatial/cytassist_image.tiff` 逐字节一致，TIFF 私有标签也包含该切片序列号和 CytAssist 标记。将 `Visium_HD_6p5mm_Human_Heart_tissue_image.btf` 用作 **显微镜 H&E 原图**，重做与发布版对应的图像配准。发布结果的低分辨率 PNG 用于显示与验证。[官方 H&E 原图下载](https://s3.us-west-2.amazonaws.com/10x.files/samples/spatial-exp/4.1.0/Visium_HD_6p5mm_Human_Heart/Visium_HD_6p5mm_Human_Heart_tissue_image.btf) · [10x 图像输入要求](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/count-visium-hd)
- 使用 `Visium_HD_6p5mm_Human_Heart_probe_set.csv`。服务器离线时，提前取得 `H1-ZX2RP3Y` 对应的 HD 载玻片布局文件，并在清单里加入 `slidefile`。如果自动配准需要人工修正，可通过 `loupe_alignment` 提供同一切片的对齐结果。[10x 载玻片与布局文件说明](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-slide-parameters)

## 在服务器上运行

把下面的清单保存为 `/data/heart/manifest.json`。文件名保持一致，路径按服务器实际位置调整。

```json
{
  "sample_id": "tenx:Visium HD:visium-hd-cytassist-6p5mm-human-heart--06669997a8ab",
  "source_sample_id": "visium-hd-cytassist-6p5mm-human-heart",
  "study_id": "tenx:Visium HD:visium-hd-cytassist-6p5mm-human-heart",
  "source_url": "https://www.10xgenomics.com/datasets/visium-hd-cytassist-6p5mm-human-heart",
  "organism": "Homo sapiens",
  "observed_sample_state": "normal",
  "platform": "Visium HD",
  "assay": "probe",
  "count_policy": "filtered_probes",
  "engine_version": "4.1.0",
  "fastqs": "inputs/Visium_HD_6p5mm_Human_Heart_fastqs",
  "fastq_sample": "Visium_HD_6p5mm_Human_Heart",
  "reference": "references/refdata-gex-GRCh38-2020-A",
  "reference_id": "GRCh38-2020-A",
  "image": "inputs/Visium_HD_6p5mm_Human_Heart_tissue_image.btf",
  "cytaimage": "inputs/Visium_HD_6p5mm_Human_Heart_image.tif",
  "probe_set": "inputs/Visium_HD_6p5mm_Human_Heart_probe_set.csv",
  "slide": "H1-ZX2RP3Y",
  "area": "A1"
}
```

从项目根目录执行。以下资源参数适用于 32 核、128 GB 的服务器。

```bash
mkdir -p /data/heart/inputs
tar -xf /data/heart/Visium_HD_6p5mm_Human_Heart_fastqs.tar -C /data/heart/inputs
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .

python -m sbasecouter recount --manifest /data/heart/manifest.json \
  --out /data/heart/recount --engine /opt/spaceranger-4.1.0/spaceranger \
  --cores 32 --memory-gb 120 --plan

python -m sbasecouter recount --manifest /data/heart/manifest.json \
  --out /data/heart/recount --engine /opt/spaceranger-4.1.0/spaceranger \
  --cores 32 --memory-gb 120 --resume
```

实际重计数的验收文件在 `/data/heart/recount`。核查 `recount.json`、`engine.log`、Space Ranger 的 `outs/web_summary.html` 和 `dataset/catalog.json`。目录中的 2、8、16 µm 三个分辨率应各自通过计数、条码、坐标和图像核查。现有 8、16 µm 发布版结果可用于对照。[项目重新计数与恢复说明](raw_recount.md) · [Space Ranger 三种分辨率的输出约定](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview)
