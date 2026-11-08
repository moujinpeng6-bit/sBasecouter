# 细胞条码、邻居与空间表达

`context` 接收一张由 `build` 或 `recount` 整理好的切片。原始实验条码继续用于匹配表达和坐标。输出另加 `spatial_observation_id`，以及能追溯输入文件、标注和参数的 `spatial_context_id`。细胞级数据保留 `spatial_cell_id`。邻居类型与表达保存在数值字段，供后续分析或模型读取。

```powershell
.\.venv\Scripts\sbasecouter.exe context --input "runs/one_slice.h5ad" --out "runs/one_slice_context.h5ad" --neighbors 8 --radius-um 50 --alpha 0.25
```

每个细胞按微米坐标寻找同切片邻居。默认最多取 8 个，限制在 50 µm 内，排除自身。有可靠标注时，可加入 `--cell-types labels.csv`，并用 `--cell-types-source` 记录标注来源。表头为 `slice_id,barcode,cell_type`，其中切片编号对应输入的 `uns["spatial_reference"]["slice_id"]`。已有 `obs["cell_type"]` 时可直接使用。缺少标签的细胞标为 `unknown`，参与距离筛选。只有已知类型相同的细胞对才触发同类型加权。

| 输出 | 含义 |
| --- | --- |
| `X` | 发布者的计数，条码及基因顺序原样保留 |
| `obsp["spatial_neighbors"]` | 半径范围内的邻居关系，可非对称 |
| `obsp["spatial_rna_weights"]` | 计算邻居 RNA 均值时实际使用的归一化权重 |
| `obsm["neighbor_type_fractions"]` | 各类型占邻居的比例，未标注类型单列，列名在 `uns["spatial_context"]["neighbor_type_columns"]` |
| `layers["neighbor_expression"]` | 邻居的每万 RNA 计数加权均值 |
| `layers["spatial_context_expression"]` | 本细胞与邻居均值的混合结果，按同一条码和基因索引保存 |

对 RNA 总计数大于零的细胞，先把 RNA 表达缩放到每万计数。有 RNA 的邻居参加表达均值。同类型邻居默认权重为 2，其他邻居权重为 1；这一选择可用 `--same-type-weight` 调整。邻居类型比例使用不加权的实际邻居数。

对本身有 RNA 也有可用邻居的细胞，默认结果为 75% 本细胞表达加 25% 邻居表达。缺少可用邻居时保留本细胞标准化结果。启用混合时，RNA 计数为零的细胞可以采用可用邻居的均值，并单独标记为仅邻居估计。`--alpha 0` 保留本细胞标准化结果。非 RNA 对照特征继续保留在原计数层。

这里借鉴 scGPT-spatial 按真实空间位置选择邻居的设计。scGPT-spatial 在模型训练中用邻域信息构建重建目标。当前实现是可复算的邻域加权基线。对于测序型平台，可在上游参考原版 scRecounter 的参考基因组、读段比对、UMI 去重和计数规则。对于 Xenium，输入是探针检测并分配给已分割细胞的计数。当前上下文步骤从发布者已有的基因 ID 和计数开始。[scGPT-spatial 论文](https://doi.org/10.1101/2025.02.05.636714) · [scGPT-spatial 源码](https://github.com/bowang-lab/scGPT-spatial) · [原版 scRecounter](https://github.com/ArcInstitute/scRecounter) · [Xenium 输出约定](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/outputs/xoa-output-overview)

| 观测 | 用于计算邻居的坐标 | 默认半径 |
| --- | --- | --- |
| Xenium 细胞 | 实验记录中的微米质心 | 50 µm |
| Visium spot | 按交错阵列行列恢复的标称微米位置，最近点中心间距 100 µm | 150 µm |
| Visium HD bin | 按原生 bin 尺寸换算的正方形网格微米位置 | bin 边长的两倍 |

Visium 和 Visium HD 的邻居限定在 `in_tissue` 的观测中，各分辨率独立运行。`obsm["context_xy_um"]` 保存用于邻域计算的位置，原有的图像像素坐标仍用于切片显示。类型比例和同类型加权只对细胞级标注启用。[10x 空间输出格式](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/spatial-outputs) · [10x 网格和点间距](https://cdn.10xgenomics.com/image/upload/v1722285090/support/documents/CG000686_VisiumHDSpatialApplicationsImagingGuidelines-TN_RevA.pdf)

## 真实切片验收

使用 10x 发布的人胰腺癌 Xenium 切片，并加入 STHELAR 发布的 `pancreatic_s0` 细胞类型分组。STHELAR 重新编号了细胞。验收按原生微米坐标一对一匹配细胞质心，逐行核对 RNA 计数与全部特征计数，再把分组写回 10x 原始条码。`Other` 归为 `unknown`。标签是研究发布的注释，保留其来源以便使用者自行判断。[10x 原数据](https://www.10xgenomics.com/datasets/pancreatic-cancer-with-xenium-human-multi-tissue-and-cancer-panel-1-standard) · [STHELAR 数据](https://doi.org/10.6019/S-BIAD2146) · [STHELAR 重新编号源码](https://github.com/MICS-Lab/STHELAR/blob/master/src/_2_annotations/_2-4_add_cells.py) · [验收统计](../examples/context_acceptance.json)
