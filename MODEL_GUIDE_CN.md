# AquaEsterNet-StructMine 模型指南

## 1. 这是什么

AquaEsterNet-StructMine 不是单一神经网络，而是一套可复现的结构挖掘和证据整合流程。它以 Bur01 结构为约束，回答：

> 哪些蛋白在整体折叠、覆盖度、催化核心和注释证据上接近 Bur01，值得进入实验验证？

项目中的三个名称各有边界：

- **AquaEsterNet**：面向水相酯形成机制和候选发现的整体框架。
- **StructMine**：负责结构搜索、链级对齐、注释、预筛分和导出的模块。
- **WASE-screen**：当前透明预筛分，不是未来正式 WASE 模型。

## 2. 输入和输出

### 输入

核心查询为 PDB 8HFW chain A：

```text
文件：data/query/Bur01m.pdb
残基：12-294
长度：283 aa
催化三联体：Ser111 / Asp241 / His274
```

目标数据可以是：

- Foldseek 数据库，例如 AlphaFold/Swiss-Prot 或 PDB。
- 本地 PDB/mmCIF 小集合。
- 由预测结构构建的自定义数据库。
- 通过 ProstT5 转换的序列集合。

### 输出

每个候选同时保留：

- 结构相似性和覆盖度。
- 查询链和目标链。
- 对齐序列和催化残基映射。
- PDB、UniProt、AlphaFold DB 或 NCBI 身份。
- reviewed 状态、pLDDT、蛋白名称、物种和序列。
- WASE-screen 分数、证据覆盖度、不确定性和缺失证据。
- 运行配置、输入哈希、软件版本、数据库指纹和错误。

## 3. 基本原理

```text
8HFW chain A
      |
      v
Foldseek 大库搜索 / TM-align 小集合对齐
      |
      v
TM-score + RMSD + 双向覆盖度过滤
      |
      v
Ser111/Asp241/His274 映射与几何检查
      |
      v
RCSB / UniProt / AlphaFold / NCBI 注释
      |
      v
WASE-screen 透明预筛分
      |
      v
运行清单 + Excel/TSV/FASTA + 网页可视化
```

### 3.1 Foldseek 搜索

Foldseek 把蛋白三维结构编码为适合快速搜索的 3Di 表征，并通过高效预过滤和结构对齐从大规模数据库中找候选。生产路线使用真实 Foldseek 输出，不用 Python 预览替代最终验收。

项目保留 Foldseek 的 `evalue`、`bits`、`prob`、TM-score、RMSD、覆盖度和序列一致性等指标，并在运行清单中记录 Foldseek 和数据库版本。

### 3.2 TM-align 预览

小型结构集合使用 `tmtools==0.3.0` 调用 TM-align。它先按结构和序列联合寻找对应关系，再计算：

- `qtmscore`：按查询长度归一化的 TM-score。
- `ttmscore`：按目标长度归一化的 TM-score。
- `rmsd`：结构叠合后的坐标偏差。
- `aligned_query` / `aligned_target`：真实对齐序列。

本地路线不生成 Foldseek 的 E-value、bit score 或 probability，因此这些字段为空。

### 3.3 链级处理

PDB 和 mmCIF 可能包含多个蛋白链。v0.3 对每条链独立处理：

- 查询链必须显式指定。
- 目标结构的每条蛋白链独立对齐。
- 不允许把不同链拼接为一条假蛋白。
- 输出 `query_chain` 和 `target_chain`。

### 3.4 催化核心映射

系统把查询的 Ser111、Asp241 和 His274 投影到目标对齐位置，记录：

- 映射后的残基编号和身份。
- 三个残基中保持同一氨基酸的数量。
- Ser-His 和 Asp-His 侧链最短距离。
- 是否落在预设的几何接触范围。

该步骤只提供催化核心证据。它不证明目标酶具有同样的底物通道、动力学或选择性。

### 3.5 注释整合

目标 ID 被解析为 PDB、UniProt 或 NCBI accession，再按配置访问：

- RCSB：实验结构、蛋白名称、物种、链和 UniProt 映射。
- UniProt：reviewed 状态、功能、基因、交叉引用和序列。
- AlphaFold DB：模型版本、pLDDT、结构和 PAE 链接。
- NCBI Protein：蛋白名称、物种和序列。

请求具有批量、重试、指数退避、磁盘缓存和错误列。最终序列来源写入 `final_sequence_source`。

### 3.6 AlphaFold2-PTM 预测层

当输入只有序列而没有结构时，v0.3 使用 LocalColabFold 1.6.1 调用 AlphaFold2-PTM。它不是把 AlphaFold 嵌入 WASE-screen 分数，而是生成一个有独立证据和清单的查询结构：

```text
单体 FASTA
  -> 格式、长度和 X 比例校验
  -> 公共 MMseqs2 MSA，按序列哈希缓存
  -> Evoformer 处理序列/MSA 表征并迭代 recycle
  -> Structure Module 生成三维坐标
  -> pLDDT、pTM 和 PAE 评估模型可信度
  -> rank_001 Amber 松弛
  -> 作为 Foldseek/TM-align 的动态查询
```

- `pLDDT` 是每残基局部可信度，不是活性分数。
- `pTM` 反映整体拓扑可信度，不证明催化构象正确。
- `PAE` 描述残基对相对位置的不确定性，可发现结构域相对取向不可靠。
- 五个标准模型共享同一输入证据但使用不同 AF2-PTM 参数；按 pLDDT 排名。

8HFW chain A 与新序列进行全局对齐，自动映射 `Ser111/Asp241/His274`。用户可以修正编号，修正值和残基身份都会进入预测清单。身份不符时不会伪造催化分。

MSA 请求和本地 GPU 推理解耦。公共请求严格串行、失败最多重试三次，默认按 15 秒、30 秒指数退避，原始服务错误进入事件和清单。GPU 通过进程锁保证一次只运行一个任务；取消信号会终止对应 WSL 进程；失效锁可恢复；显存不足不会回退 CPU。JAX 显式关闭 ColabFold 统一内存，以适配 16 GB RTX 4080；这只改变内存分配，不改变模型算法。rank_001 由 CPU OpenMM 做 Amber 松弛。

## 4. WASE-screen v0.1

### 4.1 结构适配分，45 分

```text
25 × min(qTM, tTM)
+ 10 × sqrt(qcov × tcov)
+ 10 × clamp(1 - RMSD / 5, 0, 1)
```

使用双向 TM-score 可避免长短蛋白归一化差异造成偏差；覆盖度防止只对齐一个局部片段；RMSD 项奖励更紧密的叠合。

### 4.2 催化核心分，30 分

- GXSXG-like 模体：5 分。
- 映射后的 Ser/Asp/His 身份：15 分。
- 催化几何接触：10 分。

### 4.3 注释可信度，15 分

- 酯酶、脂肪酶或水解酶相关功能注释。
- Swiss-Prot reviewed 或实验结构证据。
- AlphaFold 结构置信度。

### 4.4 实验准备度，10 分

- 可获得且不含未知 X 的蛋白序列。
- 合理的蛋白长度。
- 足够的目标覆盖度。
- 可追溯 accession。

### 4.5 层级判定

高优先级必须满足：

```text
wase_screen_score >= 70
wase_evidence_coverage >= 75
min(qTM, tTM) >= 0.50
catalytic_core_score >= 20
```

证据覆盖率不足时，即使观察到的部分分数较高，也不能进入高优先级。

### 4.6 当前缺失证据

正式 WASE 仍需要至少包括：

- 底物通道和瓶颈专一性。
- 乙醇或其他醇在活性位点的驻留。
- 水排斥和局部水化状态。
- 产物逃逸路径。
- 底物兼容性。
- 发酵生态来源先验。
- 适应性位点富集和祖先增益。

这些项目始终写入 `wase_missing_evidence`，避免把结构预筛分误写为正式机制预测。

## 5. 创新点

1. **用结构约束替代零阈值演示**：固定阳性和负对照，P00520 必须被默认阈值排除。
2. **双向 TM-score 和链级对齐**：减少长度偏差和多链拼接造成的假命中。
3. **催化残基映射进入同一记录**：不仅看整体结构，也检查 Ser/Asp/His 身份和几何。
4. **跨库证据可追溯**：注释缓存、失败原因、输入哈希和软件版本均进入运行清单。
5. **CLI 与网页共用服务层**：相同配置在两种界面下产生相同结果逻辑。
6. **明确区分 WASE-screen 与正式 WASE**：输出不确定性和缺失证据，不宣称直接预测 S/H ratio。
7. **预测与挖掘具有双向谱系**：rank_001、prediction ID、run ID、输入哈希和参数校验信息可互相追溯。
8. **动态催化参考位点**：新序列先映射 8HFW 三联体，结构搜索不再假定所有查询都使用 111/241/274 编号。
9. **隐私和 GPU 资源是显式运行状态**：MSA 传输确认、缓存命中、串行请求、GPU 锁、取消和显存错误均进入清单。

## 6. 怎么用

### 网页

```text
双击 start_ui.bat
```

网页固定使用 `http://127.0.0.1:8501`。首次双击 `install_ui_autostart.bat` 后，Windows 登录时会自动启动后台守护程序；网页进程异常退出后会自动恢复，不需要一直保留命令窗口。`ui_status.bat` 用于查看状态，`stop_ui.bat` 用于停止，`uninstall_ui_autostart.bat` 只移除登录自启动。

在左侧选择查询、链、数据库和阈值，点击“运行”。候选、结构、历史和诊断在同一工作台中查看。

“AF预测”视图支持粘贴/上传 FASTA、选择已有候选、质量档、模板开关、逐次隐私确认、催化位点预览和修正、仅预测、预测并挖掘、取消，以及 pLDDT、PAE 和三维结构检查。

在线专家演示时双击 `install_online_ui.bat`。模型和数据库仍在本机运行，外部浏览器通过带随机密码的 HTTPS 通道实时使用同一服务层。完整操作和安全限制见 `ONLINE_DEMO_GUIDE_CN.md`。

### AlphaFold2 预测

```powershell
python -m bur01_structmine.cli predict --fasta query.fasta --preset standard --confirm-public-msa

python -m bur01_structmine.cli predict --candidate-run <run_id> --accession Q9FVW3 --preset standard --confirm-public-msa

python -m bur01_structmine.cli predict-and-mine --fasta query.fasta --preset standard --confirm-public-msa --mining-config configs/alphafold_swissprot.yaml
```

标准档固定为 5 个 AF2-PTM 模型、3 次 recycle、1 个 seed、rank_001 Amber 松弛。模板只有在增加 `--templates` 时开启。

### 环境诊断

```powershell
python -m bur01_structmine.cli doctor --config configs/alphafold_swissprot.yaml
```

### 固定基准

```powershell
python -m bur01_structmine.cli benchmark --engine foldseek
```

仅做小集合预览：

```powershell
python -m bur01_structmine.cli benchmark --engine tmtools
```

### 生产搜索

```powershell
python prepare_databases.py --config configs/alphafold_swissprot.yaml --only alphafold_swissprot
python -m bur01_structmine.cli production-run --config configs/alphafold_swissprot.yaml
```

已验证参考：AlphaFold/Swiss-Prot `v6 2025_03` 数据库包含 590,183 个结构。8HFW 标准预测得到 5 个模型，rank_001 平均 pLDDT `98.07`、pTM `0.950`，相对实验结构 `min(qTM,tTM)=0.9925`；Q9FVW3 标准预测平均 pLDDT `83.54`、pTM `0.740`，相对 AlphaFold DB v6 `min(qTM,tTM)=0.7815`。预测 8HFW 的 Foldseek 基准全部通过，完整数据库搜索得到 999 个原始命中、12 个过滤候选和 12 个有效注释，零警告。

## 7. 如何解读候选

建议按以下顺序判断：

1. 检查 `min(qtmscore,ttmscore)`、RMSD 和双向覆盖度。
2. 检查催化三联体是否完整映射及几何是否合理。
3. 检查 reviewed/实验结构证据和 AlphaFold 置信度。
4. 查看 `wase_evidence_coverage` 和 `wase_missing_evidence`，不要只看总分。
5. 检查最终序列是否完整、长度是否合理、来源是否可追溯。
6. 进入表达和功能验证前，人工查看结构叠合、PAE、通道和底物位姿。

## 8. 边界和风险

- Foldseek 和 TM-align 的结构相似性不能直接推断酯化活性。
- GXSXG-like 模体可能出现在无关蛋白中，必须结合整体结构和催化核心。
- AlphaFold 的低置信度区段可能影响覆盖度和局部几何。
- 网络注释可能暂时失败，应检查 `annotation_status` 和 `annotation_errors`。
- WASE-screen 是筛选工具，不是实验结论。
