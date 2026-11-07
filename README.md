# scRecouter-spatial

把 SRAgent-spatial 验证过的**人类空间转录组**整理成可重用的 AnnData 数据集。一个输出文件对应一张切片和一种观测粒度，记录表达、坐标、图像来源和处理参数。健康、患病、癌旁、混合或未知状态沿用来源证据。

```text
SRAgent-spatial 搜索、筛选、配文件、验证对应
  → batch.json 指向同一切片的表达、坐标、比例和图像
  → scRecouter-spatial 核对文件指纹并按 ID 对齐
  → 每切片一个 .h5ad + 逐样本 catalog.json
```

## 和论文的对应

论文里 SRAgent 负责发现数据和整理元数据。scRecounter 接收原始测序读段，用 STARsolo 对齐参考基因组并统一计数。scBaseCount 把整理后的数据与来源信息发布成仓库。本项目 v0.1 从发布者已经处理好的计数起步，v0.2 给细胞级数据增加空间邻域表达。v0.3 新增 `recount` 和 Nextflow 批量入口，调用空间平台的计数工具，再复用空间对应检查和 AnnData 输出。[原始数据重处理](docs/raw_recount.md) · [Cell 论文](https://doi.org/10.1016/j.cell.2026.08.025) · [原版 scRecounter](https://github.com/ArcInstitute/scRecounter)

**标准化和归一化各管一件事。** 标准化统一文件结构、条码对齐、坐标单位、质量统计和来源记录。归一化是可选的数值变换。基础版一直把发布者提供的计数保存到 `X`，方便复查和重新选择分析方法。

## 使用

Python 3.11 或更新版本。Windows PowerShell 示例。

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\sbasecouter.exe build --report "E:/path/to/batch.json" --out "runs/human_spatial"
```

- `--report` 可重复，输入 SRAgent-spatial 的 `research-spatial` 所产生的本地 `batch.json`。`--sample-id` 可用报告中的精确 ID 选切片。
- `--normalize-total 10000` 增加一个 `log1p_total` 层。每个观测按 RNA 基因计数之和缩放到指定值，然后取 `log(1 + x)`。零基因计数仍为零。对照探针和其他非 RNA 特征保留在 `X`，归一化层中相应列为零。不同平台和不同基因面板的结果分别使用。
- 输出目录须为新目录。逐样本写入结果，失败会记录原因并继续处理其他样本。原图留在源目录，`catalog.json` 记录图像路径、来源链接和 SHA-256。跨机器使用时需带上原图并重新指向其位置。

## 文件里有什么

| 字段 | 内容 |
| --- | --- |
| `X` | 发布者提供的稀疏整数计数，行是观测，列是特征 |
| `obs` | 原始条码、总计数、RNA 基因计数和检测到的特征数。Visium 还保留组织内标记与阵列行列 |
| `var` | 发布者特征 ID、名称、类型和可用的基因组字段 |
| `obsm["spatial"]` | 切片的原生 X/Y 坐标 |
| `obsm["image_xy"]` | 与该切片指定图像对应的 X/Y 像素坐标 |
| `uns` | 空间单位、图像尺寸、坐标比例、来源和输入文件哈希 |

| 平台 | 一行代表 | 原生坐标 |
| --- | --- | --- |
| Visium | 一个 spot | 全分辨率图像像素，按来源比例映射到低分辨率或高分辨率 PNG |
| Visium HD | 一个 bin | 全分辨率图像像素，保留原生 bin 尺寸，要求完整计数和对应的 Parquet 坐标 |
| Xenium | 一个已分割细胞 | 微米，按实验文件的每像素微米数映射到形态图像 |

只有通过来源空间检查并带有完整表达的切片才进入正式输出。只有矩阵头的样本会以 `full_count_matrix_required` 记录。每张切片使用自己的参考系，跨切片的坐标各自保留。[Visium 坐标约定](https://www.10xgenomics.com/support/cn/software/space-ranger/2.1/tutorials/outputs/spatial-outputs) · [Xenium 坐标约定](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html) · [AnnData 格式](https://anndata.readthedocs.io/en/stable/fileformat-prose.html)

## 验收

v0.2 可对细胞级切片计算邻居图、细胞类型比例、邻居表达和空间上下文表达。人胰腺癌 Xenium 的类型标注与原条码逐个核对后已实测。命令及计算方式见 [空间上下文](docs/spatial_context.md)，核对结果见 [v0.2 验收统计](examples/context_acceptance.json)。

| 已重读的真实切片 | 数据量 | 验收范围 |
| --- | --- | --- |
| 人肺 Visium，两张切片 | 4,034 和 4,882 个 spot | 重新打开结果，逐项比对源计数、条码、基因 ID，并抽查归一化 |
| 人胰腺癌 Xenium，一张切片 | 190,965 个细胞、474 个 RNA 基因特征 | 同样核对完整计数，保留来源中的非 RNA 特征和细胞坐标 |
| 健康人心脏 Visium HD，一张切片 | 10,615,406 个原生 2 µm bin | 复用同一官方归档的完整矩阵及配套特征表，重新读取输出并核对尺寸、计数条目和空间比例 |

心脏原搜索报告的矩阵头被正确跳过，补入同一发布归档的完整矩阵后才生成结果。验收记录只包含公开来源与统计。[验收详情](examples/acceptance.json)

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```
