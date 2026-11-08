# 从原始数据重新生成空间表达矩阵

原版 scRecounter 先用部分读段试选参数，再用完整读段比对和计数。空间版沿用「明确输入、固定计数规则、保存矩阵和来源」的流程。平台与实验方法由样本清单明确指定。清单中的 FASTQ、切片图像、捕获区及参考库需要预先配齐。[原版 scRecounter](https://github.com/ArcInstitute/scRecounter) · [STARsolo 官方说明](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)

| 平台和方法 | 重新处理工具 | 可选择的计数规则 | 矩阵行 |
| --- | --- | --- | --- |
| Visium 或 Visium HD，3′ 测序 | Space Ranger | `exonic` 排除内含子读段，`exonic_intronic` 包含内含子读段 | spot 或 bin |
| Visium 或 Visium HD，探针法 | Space Ranger | `filtered_probes` 按探针表过滤，`all_nondeprecated_probes` 纳入未弃用探针 | spot 或 bin |
| Xenium，原位检测 | Xenium Ranger | `resegment` 重做细胞分割并分配已定位转录本，可指定 `expansion_um` | 重新分割的细胞 |

Space Ranger 的两种 3′ 规则与原版 `Gene` / `GeneFull` 对比的是类似的生物学取舍，具体比对、分配和去重由各自的程序执行。探针法的两个选项改变纳入的探针集合。Visium HD 分别验证 2、8、16 µm，保留每个分辨率独立的表达条码、位置、图像和比例因子。[Space Ranger 参数](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/command-line-arguments) · [Space Ranger 输出](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview) · [Xenium Ranger 重新分割](https://www.10xgenomics.com/support/software/xenium-ranger/latest/tutorials/running-pipelines/XR-resegment)

## 准备一张人类切片

保存一个 JSON 清单。以下内容演示普通 3′ Visium，需替换成实际样本。路径可相对于该清单，也可以使用绝对路径。FASTQ 目录只放指定样本，所有 lane 和 chunk 的 R1、R2 都要配齐。实际分析由指定版本的 Space Ranger 读取整组 FASTQ。[Space Ranger FASTQ 命名](https://www.10xgenomics.com/support/software/space-ranger/latest/tutorials/inputs/fastqs-specifying-fastqs)

```json
{
  "sample_id": "exact-sample-id-in-sragent-report",
  "source_sample_id": "published-sample-id",
  "study_id": "published-study-id",
  "source_url": "https://example.org/replace-with-the-published-sample-page",
  "source_report": "../discovery/batch.json",
  "organism": "Homo sapiens",
  "platform": "Visium",
  "assay": "three_prime",
  "count_policy": "exonic_intronic",
  "engine_version": "4.1.0",
  "fastqs": "../inputs/fastq",
  "fastq_sample": "HumanSlice",
  "reference": "../reference/refdata-gex-GRCh38-2024-A",
  "reference_id": "GRCh38-2024-A",
  "image": "../inputs/tissue-image.tif",
  "slide": "replace-with-slide-serial",
  "area": "A1"
}
```

- `source_report` 可选。提供时按 `sample_id` 查找唯一记录，核对来源样本号、研究号、物种与平台，并带入报告中的疾病状态和来源证据。清单仍需给出与该切片对应的原始输入路径。
- Visium HD 将 `platform` 改成 `Visium HD`，并提供 `cytaimage`。探针法将 `assay` 改成 `probe`，同时提供 `cytaimage`、`probe_set` 并选择探针计数规则。离线分析可通过 `slidefile` 指向匹配的官方切片布局文件。需要手动图像配准时，可提供同一切片的 `loupe_alignment` 文件。
- Xenium 将 `platform`、`assay`、`count_policy` 改成 `Xenium`、`in_situ`、`resegment`，用 `xenium_bundle` 指向完整 XOA 输出目录。`engine_version` 改成实际安装的 Xenium Ranger 版本，例如 `4.0.0`。重新分割后应重新核查细胞类型标签与新细胞 ID 的对应。
- 比较同一切片的两种计数规则时，用两个清单文件，保持切片 ID 和源文件相同，分别指定 `count_policy` 并写到不同输出目录。

## 运行与读取结果

实际重新计数需要装有对应 10x 工具的 Linux 环境。批量模式需要 Nextflow 24.10 或以上版本。下面的命令从项目根目录运行。Windows 也可先用 `--plan` 检查清单与路径。

```bash
python -m pip install -e .
python -m sbasecouter recount --manifest manifests/slice-a.json --out runs/slice-a --plan
python -m sbasecouter recount --manifest manifests/slice-a.json --out runs/slice-a --cores 16 --memory-gb 128
nextflow run workflows/recount.nf --manifests manifests --outdir runs/batch --cores 16 --memory_gb 128
```

单样本 `recount` 默认使用新目录。v0.4 可用 `--resume` 重开同一个输出目录，恢复前核对原始文件、参考库、引擎版本和运行参数。厂商引擎中断时，用相同命令继续原有运行。引擎已完成时，继续空间核查并复用已校验的矩阵。Nextflow 批量入口自动传入 `--resume`，其自身缓存关闭，由项目的校验记录判断可复用项。单样本失败后，Nextflow 仍会调度其他样本并在结束时返回失败状态。引擎仍在运行或锁文件待处理时，请根据厂商日志确认运行状态。[Nextflow 执行与资源](https://training.nextflow.io/latest/config_exec/02_resources_and_retries/) · [Space Ranger 恢复说明](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/sr-troubleshooting)

任务先写 `recount.json` 记录状态。进入计数时还会写 `inputs.json` 和 `engine.log`。生成的矩阵必须先通过条码、坐标、图像和尺度核对，才进入 `dataset/catalog.json` 的成功记录。`.h5ad` 的 `X` 保存重新处理的整数计数，`source_record_json` 保存引擎版本、规则、参考库标识和来源关联。对细胞级 Xenium 输出，可接着运行 [`context`](spatial_context.md) 生成带邻居信息的另一层表达矩阵。

当前验收覆盖命令参数、来源编号匹配、漏读段时的预检、模拟 10x 输出到 AnnData 的转换，以及 HD 缺一个分辨率时的记录。真实原始数据运行需要 Linux 或 WSL、对应版本的 10x 引擎和同一切片的完整输入。
