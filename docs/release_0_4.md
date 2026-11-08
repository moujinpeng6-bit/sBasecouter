# v0.4 运行入口和验收范围

## 从 SRAgent-spatial 到表达矩阵

1. `build` 从 SRAgent-spatial 的批量报告整理已核查的发布者矩阵。`--resume` 校验并保留已有成功条目，可以追加新的报告。再次运行时需包含之前的样本记录。
2. `prepare` 接收发现报告与原始文件清单，按研究、样本、平台和人类物种证据精确匹配。支持公开文件断点续传，并按大小和 SHA-256 验证。来源记录无法定位原始文件时，在 `preparation.json` 中标为 `raw_recipe_required`。
3. `recount` 或 Nextflow 在具备厂商引擎的 Linux 上处理准备好的输入。`--resume` 对比输入和命令。重计数完成后核查条码、坐标、尺度和图像，再写 `.h5ad`。成功样本在下次运行时按 SHA-256 复用。
4. `context` 读取同一切片的 `.h5ad`。细胞、spot 和 bin 分别构建邻域，并将衍生表达保存为独立层。

```bash
python -m sbasecouter prepare --report batch.json --recipes raw-recipes.json --out prepared
nextflow run workflows/recount.nf --manifests prepared/manifests --outdir recounted
python -m sbasecouter build --report original/batch.json --report added/batch.json --out catalog --resume
```

`raw-recipes.json` 的顶层为 `samples` 数组。每条含 `sample_id`、`source_sample_id`、`study_id` 和 `platform`，取值与报告一致。`manifest` 保存计数所需字段，格式见[原始数据重处理](raw_recount.md)。`files` 的键是样本目录中的目标相对路径，例如 `fastqs/HumanSlice_S1_L001_R1_001.fastq.gz`。每个值含直链 `url`、整数 `bytes` 和小写 `sha256`。`reference` 相对于原始文件清单的位置解析，指向已解包的参考库。下载归档需按来源说明解包并核对成员后，再生成单样本清单。

## 验收边界

六张真实人类切片使用已发布计数完成空间邻域计算，覆盖肺、淋巴结、胰腺和心脏。心脏为完整 2 µm 原生网格。恢复测试在第一张真实切片建好后中断进程，续跑得到全部三张切片，第一张矩阵的哈希与修改时间一致。厂商重计数的命令调度、恢复和 2、8、16 µm 输出检查使用模拟引擎测试。[逐样本验收记录](../examples/release_0_4_acceptance.json)

验证真实重计数还需要受支持的 Linux 环境、足够的内存和磁盘、厂商工具以及同切片的完整原始文件。对于 `resegment`，需要完整 Xenium 原始分析包。Visium 和 Visium HD 需要完整 FASTQ、参考库和与切片相符的图像及布局资料。[Space Ranger 系统要求](https://www.10xgenomics.com/support/software/space-ranger/latest/release-notes/space-ranger-system-requirements) · [Xenium Ranger 系统要求](https://www.10xgenomics.com/support/software/xenium-ranger/latest/release-notes/XR-system-requirements)

跨队列验证需要独立的人类样本集和明确的健康或患病证据。邻域表达是可复算的估计，保留原始整数计数供复核。
