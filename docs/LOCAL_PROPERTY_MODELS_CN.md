# AquaEsterNet 本地性质预测

输入一个或多个 AquaDesign 批次，以及分别命名的底物 SMILES。程序从各批次亲本坐标提取序列，核对哈希值、候选标识、突变位置和替换数。相同序列只推理一次，每个候选仍与产生它的亲本比较。模型在本机运行，结果保存在 results/properties 下。

| 模型 | 输入与输出 | 解释 |
|---|---|---|
| SWI / SoDoPE | 氨基酸序列 → 组成加权指数、经原作者 logistic 函数转换的 0–1 分数 | 表示大肠杆菌可溶表达倾向。标签、载体、培养及诱导条件需在表达实验中控制。[原作者实现](https://github.com/Gardner-BinfLab/SoDoPE_paper_2020/tree/c28f26d94aa6285ff08270709fc8a95806ddb1db/SWI) |
| TmProt 1.0 | ESM-2 650M 加作者发布的 LoRA 和回归头 → 预测熔解温度（°C） | 运行时核对适配器与回归头权重。上限按作者 512 个 token 设置。超过 510 个残基的序列将 Tm 留空，并注明原因。近缘突变间的小数差异需实验鉴别。[原作者代码](https://github.com/loschmidt/TmProt/tree/02d4c12b39ffae8a6cde64a401210390ea9cb23f/tmprot-1.0) · [ESM-2 基础模型](https://huggingface.co/facebook/esm2_t33_650M_UR50D) |
| EnzRank | 序列和一种小分子的含手性 Morgan 指纹（半径 2，2048 位） → 酶–底物配对分数（0–1） | 使用作者发布的最终 CNN、原氨基酸编码和左填充。辛酸与乙醇分别输入。模型未联合描述两者在水相中的反应，也未估计酯合成速率、水解速率或 S/H。[原作者代码与权重](https://github.com/maranasgroup/EnzRank/tree/bf9d5388e07d719b94bad3aabc010b78bcd1e2c5) |

## 安装与运行

推荐使用项目的 Python 3.11 环境执行。首次安装需要约 2.7 GB 模型文件，另需独立 Python 运行环境和依赖空间。安装操作下载固定版本，按 SHA-256 核对后才替换本地文件。正常推理从本地文件读取，向模型传入序列与 SMILES 的过程无需访问外部服务。

~~~sh
python -m bur01_structmine.cli properties-setup --install-runtime
python -m bur01_structmine.cli properties --run results/designs/my_design --substrate "octanoic_acid=CCCCCCCC(=O)O" --substrate "ethanol=CCO"
python -m bur01_structmine.cli properties-web
~~~

重复提供 --run 可以合并批次。省略 --substrate 时使用上述两种分子。显式提供 --output 时需选择尚不存在的目录。输出中保留原输入以及 RDKit 规范化的 SMILES，不自动改变质子化状态。可用 AQUA_PROPERTY_MODELS 指定权重目录，用 AQUA_PROPERTY_PYTHON 指向已安装 requirements-properties.txt 的解释器。默认 TmProt 和 EnzRank 均使用 CPU；若已配好 CUDA 版 PyTorch，可用 --device cuda 指定 TmProt 使用显卡。

## 网页和原始结果

Windows 双击 start_properties_web.bat，或使用 properties-web 启动仅监听本机的网页。在“本地计算设置”选择批次和底物，点击重新计算即可生成新的记录。在完整工作台中，也可以从“结构 → 本地酶序列设计”打开同一套性质预测。结果页支持切换指标和底物、筛选替换数、按候选查看原始分数，也可下载离线网页用于答辩演示。

- protein_predictions.tsv 每条唯一序列一行，包含亲本身份、序列、SHA-256、SWI 指数和可溶表达分数、Tm 及其可评估状态。
- substrate_predictions.tsv 每个序列–底物对一行，保留输入 SMILES、规范 SMILES 及 EnzRank 原始分数。
- candidate_comparison.tsv 每个设计批次–候选–底物组合一行。差值按“候选减去对应亲本”计算。亲本序列与 8HFW A 链相同时，按结构记录逐个映射成共同的 8HFW 编号，并保留原始坐标与输入序列编号。其他亲本使用自身的输入序列编号。
- manifest.json 记录输入和结果哈希、模型源版本、权重指纹、执行设备、运行软件版本、计算数量及底物输入。
- index.html 包含本批真实模型结果及图表逻辑。可在无网络的浏览器中交互查看；重新推理使用上面的本地工作台。

溶解表达分数、预测 Tm 与底物配对分数的单位和训练任务不同。界面按用户选定的单项指标排序，保留每种底物的独立结果。
