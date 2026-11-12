# Windows 原始输入与空间表达验收

## 心脏原始输入

健康人心脏 Visium HD 的显微镜 H&E 原图和 `GRCh38-2020-A` 参考库位于 `E:/scRecouter-spatial/v0.4/downloads`。这个目录也包含已核对的 FASTQ、CytAssist 图像和探针表。[下载与内容核查记录](../examples/windows_heart_inputs.json) · [Space Ranger 运行所需输入](heart_linux_recount.md) · [10x 样本发布页](https://www.10xgenomics.com/datasets/visium-hd-cytassist-6p5mm-human-heart)

下面两条命令重读已有文件并核对字节数、官方分段 ETag、SHA-256 和内部结构。缺少文件时按段取得。参考库以压缩包保存，Linux 运行前需要解包。

```powershell
.\.venv\Scripts\python.exe examples/prepare_heart_inputs.py --out E:/scRecouter-spatial/v0.4/downloads --role he
.\.venv\Scripts\python.exe examples/prepare_heart_inputs.py --out E:/scRecouter-spatial/v0.4/downloads --role reference
```

## 搜索后生成原始输入清单

脚本调用 SRAgent-spatial 的公开网页索引，从组织和平台搜索 10x 数据卡，再查找发布者的 Space Ranger 运行摘要。程序核对载玻片、捕获区、FASTQ 样本前缀及 CytAssist 图像内的载玻片标签。清单记录每项原始输入的官方链接、字节数和 ETag，并记录实际抽读的 FASTQ 归档条目。可用摘要和文件取决于发布者提供的内容。

```powershell
..\SRAgent-spatial\.venv-agent\Scripts\python.exe examples/discover_raw_visium.py --sragent-source ../SRAgent-spatial --tissue breast --keywords "Fixed Frozen" --search-label "Human Breast Cancer" --out raw_inputs.json
```

实测找到此前未验收的人类乳腺癌切片。输入只包含组织和样本描述。程序查到载玻片 `H1-BYFBQ82`、捕获区 `A1` 和 Space Ranger `3.1.1`，并登记 FASTQ 归档、两种图像、探针表和参考库。状态 `public_sources_verified` 表示远端来源和抽查的身份字段通过。归档只做了小范围读取。运行 Space Ranger 前，需要取得完整文件并准备本地路径。[逐项清单](../examples/windows_raw_discovery.json) · [10x 乳腺癌数据卡](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-breast-cancer-fixed-frozen)

## 邻居表达对皮层分层识别的影响

使用公开人类背外侧前额叶皮层 Visium 数据。14 位供体各一张切片，7 位为研究对照，7 位有精神分裂症。作者对其中 20,677 个 spot 做了组织学分层标注，全部与表达条码、坐标及图像配套文件核对。任务是根据表达预测 spot 的皮层分层或白质标注。[公开数据](https://zenodo.org/records/6921620) · [作者的供体分组](https://github.com/khodosevichlab/Schizophrenia20/blob/main/Notebook6_Visium.md#L537-L538)

每轮留一位供体测试，其他供体用于选取 512 个变异基因和训练同一种分类器。两组特征分别使用自身表达和 `75% 自身表达 + 25% 邻居平均表达`。邻居来自同一切片内 150 µm 内最多 8 个有 RNA 计数的 spot。先按 spot 的总基因计数归一化，再对两组特征取 `log1p`。宏平均 F1 给测试切片中出现的各个分层同等权重。表中是各供体分数的平均值。

| 供体分组 | 自身表达准确率 | 加邻居准确率 | 自身表达宏平均 F1 | 加邻居宏平均 F1 | F1 提升的供体 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 全部 | 28.6% | 37.4% | 0.283 | 0.335 | 12/14 |
| 对照 | 27.1% | 36.1% | 0.316 | 0.368 | 5/7 |
| 精神分裂症 | 30.0% | 38.6% | 0.250 | 0.301 | 7/7 |

这组结果衡量同一公开队列的分层识别。其他器官、平台和疾病判别需要各自验证。[逐供体分数](../examples/windows_layer_benchmark.json) · [输入文件指纹](../examples/windows_layer_inputs.json)

```powershell
.\.venv\Scripts\python.exe examples/prepare_layer_benchmark.py --sragent-source ../SRAgent-spatial --out runs/windows_benchmark
.\.venv\Scripts\python.exe examples/benchmark_spatial_layers.py runs/windows_benchmark/input.json runs/windows_benchmark/result.json
```
