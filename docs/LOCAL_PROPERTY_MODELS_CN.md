# AquaEsterNet 本地性质预测

可提交自己的蛋白序列和底物 SMILES，也可选取一个或多个 AquaDesign 批次。设计批次从亲本坐标提取序列，逐个核对序列、突变位置和替换数。同批相同序列只计算一次。结果保存在 results/properties 下。

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
python -m bur01_structmine.cli properties --fasta my_sequences.fasta --substrate "ethyl_acetate=CCOC(=O)C"
python -m bur01_structmine.cli web
~~~

重复提供 --run 可以合并批次。省略 --substrate 时使用上述两种分子。显式提供 --output 时需选择尚不存在的目录。输出中保留原输入以及 RDKit 规范化的 SMILES，不自动改变质子化状态。可用 AQUA_PROPERTY_MODELS 指定权重目录，用 AQUA_PROPERTY_PYTHON 指向已安装 requirements-properties.txt 的解释器。默认 TmProt 和 EnzRank 均使用 CPU；若已配好 CUDA 版 PyTorch，可用 --device cuda 指定 TmProt 使用显卡。

## 网页和原始结果

Windows 双击 **start_ui.bat** 或 **start_model_web.bat**，命令行可运行 **web**。默认进入 http://127.0.0.1:8501 的 StructureMine“运行”页。“序列设计”和“性质预测”位于“结构”标签旁，使用同一个侧栏和所选挖掘批次。启用了访问密码的工作台使用原密码登录。旧的独立结果链接可用 **properties-web** 在 8765 端口打开。首次安装结构设计模型的步骤见[本地酶设计说明](NATIVE_ENZYME_DESIGN_CN.md)。

1. **性质预测** 接收粘贴的单条序列、多条 FASTA 或 UTF-8 FASTA 文件。上传文件优先于文本框。每条记录使用独立名称，允许 1–2500 个标准氨基酸，一批最多 100 条。默认按序列给出三项原始预测。若输入亲本及其等长替换体，可勾选“以首条序列为亲本”，或在命令行加 --first-as-parent。FASTA 替换编号从首条序列第 1 位起算。
2. **序列设计** 可选择 8HFW、侧栏当前批次中有本地坐标和完整序列的亲本，或上传 PDB、CIF、mmCIF 并指定链。生成候选后点击“用本批候选预测性质”。“设计批次”也能合并已经生成的批次，分别采用各自的结构亲本。装有原两批 Bur01 设计的电脑可一键选择其 24 条候选。
3. 每行输入一种名称和 SMILES，点击 **提交本地预测**。网页按提交顺序在本机运行任务，关闭浏览器标签后服务进程继续计算。**任务与结果** 自动更新进度，结果地址包含任务编号，可刷新或收藏。重新启动服务后，已保存的报告可以继续打开；排队和计算中的任务需要重新提交。
4. **候选** 页可从当前筛选结果选一条记录，点击“送入性质预测”带入其注释蛋白序列。**NR挖掘** 使用所选命中的目标序列。交接时会显示输入长度；有疑问可核对 FASTA 后再提交。例如，8HFW 的数据库条目可以指向较长的 UniProt 序列，设计批次用的亲本则取自有坐标的链段。

- protein_predictions.tsv 每条唯一序列一行，包含亲本身份、序列、SHA-256、SWI 指数和可溶表达分数、Tm 及其可评估状态。
- substrate_predictions.tsv 每个序列–底物对一行，保留输入 SMILES、规范 SMILES 及 EnzRank 原始分数。
- 存在亲本对照时，candidate_comparison.tsv 每个候选与底物组合一行，差值按“候选减对应亲本”计算。设计批次的亲本与 8HFW A 链相同时，另按结构记录映射到共同的 8HFW 编号。FASTA 任务保存规范化的 input.fasta 以及每条输入名称到序列摘要的对应关系，含有相同序列的不同记录也能分别查看。
- manifest.json 记录输入和结果哈希、模型源版本、权重指纹、执行设备、运行软件版本、计算数量及底物输入。
- index.html 包含本批真实模型结果及图表逻辑。可在无网络的浏览器中交互查看；重新推理使用上面的本地工作台。

溶解表达分数、预测 Tm 与底物配对分数的单位和训练任务不同。界面按用户选定的单项指标排序，保留每种底物的独立结果。
