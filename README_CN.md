# AquaEsterNet-StructMine v0.4 中文说明

**AquaEsterNet-StructMine** 是一个面向 Bur01-like 酶的结构与序列联合挖掘工具。它以 [PDB 8HFW](https://www.rcsb.org/structure/8HFW) chain A 为查询结构，用 Foldseek/TM-align 发现结构候选，用 NCBI NR 和 MMseqs2 扩增天然序列空间，再通过本地 AlphaFold2、结构复核、跨库注释、WASE-screen 和可追溯清单形成分层证据链。

本项目不是一个直接预测酯化选择性或 S/H ratio 的黑箱模型。它解决的是更前置的问题：从结构库中找出值得进行表达、酶学和底物验证的 Bur01-like 候选。

## 最快开始

在工作台的“结构”页打开“本地酶序列设计”，可围绕 Bur01 类亲本设计少量替换，比较配体原子对序列偏好的影响，并将候选送入 AF 预测或导入结构复核。安装、算法和结果解释见[本地酶序列设计说明](docs/NATIVE_ENZYME_DESIGN_CN.md)。

双击项目根目录中的：

```text
start_ui.bat
```

它会启动后台守护程序并打开固定地址 `http://127.0.0.1:8501`。网页不再依赖命令窗口持续打开；如果 Streamlit 意外退出，守护程序会自动重启。默认界面为中文，可切换英文。

首次启用登录自启动时，双击：

```text
install_ui_autostart.bat
```

安装后，Windows 用户登录时会在后台启动网页，不会自动暴露到局域网或公网。其他控制入口：`stop_ui.bat` 停止网页，`ui_status.bat` 查看状态，`uninstall_ui_autostart.bat` 只移除自启动任务而不删除项目或结果。

命令行入口：

```powershell
# 环境诊断
python start_here.py --action doctor

# 固定基准的真实 Foldseek 验收
python start_here.py --action benchmark

# 完整 AlphaFold/Swiss-Prot 生产搜索
python start_here.py --action production-run

# 运行隔离测试
python start_here.py --action validate

# 检查项目文件是否完整
python start_here.py --action check-files

# 单体 AF2 预测。首次或缓存未命中时必须确认公共 MSA 传输
python -m bur01_structmine.cli predict --fasta query.fasta --preset standard --confirm-public-msa

# 预测 rank_001 后直接作为完整 Swiss-Prot 挖掘查询
python -m bur01_structmine.cli predict-and-mine --fasta query.fasta --preset standard --confirm-public-msa --mining-config configs/alphafold_swissprot.yaml

# 真实 NCBI NR 命中子集；需要联系邮箱和本次明确确认
python -m bur01_structmine.cli nr-demo-build --query-fasta data/query/Bur01m_seed.fasta --email name@example.org --confirm-ncbi

# 从一次结构运行扩增 NR 序列
python -m bur01_structmine.cli nr-expand --source-run <run_id> --config configs/nr_local_demo.yaml

# 导出完整 NR 服务器/HPC 作业包
python -m bur01_structmine.cli nr-export-job --source-run <run_id> --config configs/nr_full_server.yaml

# 导回并逐文件校验完整 NR 结果
python -m bur01_structmine.cli nr-import-job --bundle <nr_run_id>_server_result.zip

# NR 最近邻新颖性审计
python -m bur01_structmine.cli nr-novelty --fasta designed.fasta --backend ncbi_remote --email name@example.org --confirm-ncbi
```

`target-run` 仍保留兼容，但现在与 `benchmark` 相同，运行真实 Foldseek 固定基准，不再把 P00520 当作阳性结果。

## 它有什么用

1. 从 AlphaFold/Swiss-Prot 或 PDB 中发现结构上接近 Bur01 的候选蛋白。
2. 在序列一致性不高时，利用整体折叠、覆盖度和催化核心保守性保留潜在远缘同源物。
3. 把 Foldseek、TM-align、UniProt、AlphaFold DB、RCSB 和 NCBI 信息整理成同一张候选表。
4. 输出可用于后续基因合成、表达纯化和酶学实验的 Excel、TSV 和 FASTA 文件。
5. 为未来正式 WASE 模型准备透明、可审计的结构和证据字段。

## NR 序列扩增与新颖性审计

v0.4 参考 VenusMine 的公开 MMseqs2 搜索与聚类步骤，固定采用两轮迭代、`max-seqs=2500`、`search-type=1` 和 50% 序列一致性聚类。Foldseek 结果现在保存目标序列；种子按 accession 和序列 SHA-256 去重，聚类代表序列再直接对齐 8HFW chain A，分别映射 Ser111、Asp241 和 His274。

本机只保存一次真实 NCBI 远程 `nr` 搜索产生的命中子集，并明确标记为 `nr_remote_snapshot`。完整 NR 通过无凭据 Bash/Slurm/PBS 作业包在服务器或 HPC 运行，普通搜索不会自动下载 NR。导回结果必须通过范围、版本、路径安全和 SHA-256 校验。

NR 序列证据与结构证据严格分开。任何尚未完成 AF2 和 Foldseek/TM-align 复核的候选均显示：

```text
wase_tier = Not evaluated - structure required
wase_scope = not_evaluated_sequence_only
```

“NR挖掘”网页视图支持本机子集、服务器作业导出/导入、证据筛选和送入 AF2 复核。详细说明见 [`docs/NR_WORKFLOW_CN.md`](docs/NR_WORKFLOW_CN.md)。`configs/ncbi_prostt5.yaml` 仅保留为小型自定义 FASTA 的 legacy 入口，不代表完整 NR。

## AlphaFold2 预测功能

v0.3 起在现有 `Ubuntu-Foldseek` WSL2 中加入 LocalColabFold 1.6.1 和 AlphaFold2-PTM；v0.4 进一步允许从 NR 运行和 accession 直接选择代表序列，形成 NR → AF2 → Foldseek/TM-align 的双向可追溯谱系。

### 支持的输入

- 单条 FASTA，长度 20-1000 aa。
- 只允许 20 种标准氨基酸和 `X`；`X` 不得超过 5%。
- 可直接选择某次挖掘运行中的已有候选 accession。
- 当前 AF2 功能只支持单体蛋白，不支持复合物、配体、核酸或共价修饰。

### 标准参数

| 质量档 | 模型 | recycle | seed | Amber 松弛 |
|---|---:|---:|---:|---:|
| 快速 | 1 | 3 | 1 | 0 |
| 标准，默认 | 5 | 3 | 1 | rank_001 |
| 高采样 | 5 | 12 | 3 | rank_001 |

模型类型固定为 `alphafold2_ptm`。模板默认关闭；只有用户显式开启时才使用公共模板搜索。预测不会静默回退 CPU，JAX 未检测到 GPU 或显存不足时会明确失败。

### 两阶段运行和隐私

1. 公共 ColabFold MMseqs2 服务生成 MSA。
2. RTX 4080 在本地运行 AlphaFold2-PTM 推理；模型排名后由 CPU OpenMM 对 rank_001 做 Amber 松弛，避免与仍驻留的模型权重争用显存。

MSA 按标准化序列 SHA-256 缓存。所有公共 MSA 请求严格串行，缓存未命中时必须逐次确认序列将发送到第三方服务。网络或服务端失败最多重试 3 次，默认按 15 秒、30 秒指数退避，并把原始服务错误写入事件和失败清单。项目不下载约 940 GB 的本地 ColabFold 数据库。

在 `Ubuntu-Foldseek` 中、从挂载的项目根目录安装固定运行环境和权重：

```bash
sudo bash scripts/install_colabfold_wsl.sh
sudo bash scripts/install_alphafold_params.sh
```

权重归档和五个 PTM 模型文件均记录 SHA-256。AlphaFold 参数使用其随附的 CC BY-NC-SA 4.0 许可，不会进入 Git 仓库或发布包。

### 预测结果

每个任务写入 `results/predictions/<prediction_id>/`，包括：

- 标准化 FASTA、序列哈希和 A3M。
- 5 个标准档模型、排名、rank_001 PDB 和模型参数校验值。
- 每残基 pLDDT、pTM、PAE、日志和进度事件。
- 运行时版本、参数、MSA 缓存状态、警告和失败原因。
- 8HFW `Ser111/Asp241/His274` 的全局序列映射及人工修正记录。
- 与后续 Foldseek 运行双向关联的 prediction ID 和 run ID。

如果催化位点无法映射或残基身份不符，预测和结构搜索仍继续，但催化核心证据保持为空并明确记录缺失项。

## 查询结构

| 项目 | 内容 |
|---|---|
| PDB | 8HFW |
| 查询链 | A |
| 文件 | `data/query/Bur01m.pdb` |
| 残基范围 | 12-294 |
| 链长度 | 283 aa |
| 催化残基 | Ser111 / Asp241 / His274 |
| SHA-256 | `9f51637496713b18d3ec8f0586a1aef19d3c8f96bd3e48fda9e9af6e88da24fe` |

每次运行都会把这些信息写入 `run_manifest.json`。多链结构不会再被拼接成一条假结构，每条链独立比较并输出 `query_chain` 和 `target_chain`。

## 固定基准集

| 类型 | 对照 | 用途 |
|---|---|---|
| 自比对阳性 | 8HFW chain A | 检查 TM-score 和 RMSD 上限 |
| 相关阳性 | 1XKL、6TJ2、Q6RYA0、Q9FVW3 | 检查 Bur01-like 结构排序 |
| 长度匹配负对照 | P60174 | 排除仅由长度相近造成的假阳性 |
| 远缘负对照 | P00520 | 排除旧版错误阳性 |

验收条件：

- 8HFW 自比对 TM-score 大于 0.98，RMSD 接近 0。
- 五个阳性整体排在两个负对照之前。
- P00520 不得通过默认结构阈值。
- 最终正式验收必须由 Foldseek 完成，TM-align 预览不能替代。

基准来源和下载校验值记录在：

```text
data/benchmark/benchmark_manifest.yaml
data/benchmark/download_manifest.json
```

## 两种结构搜索引擎

### Foldseek

用于固定基准的正式验收和 AlphaFold/Swiss-Prot 大库搜索。默认 WSL 配置：

```text
发行版：Ubuntu-Foldseek
Windows 安装位置：可配置的 WSL 导入目录
Foldseek：/opt/foldseek/bin/foldseek
数据库：<项目目录>/databases
```

完整环境安装脚本：

```powershell
powershell -ExecutionPolicy Bypass -File scripts\setup_wsl_foldseek.ps1
```

完整验收脚本：

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_production_acceptance.ps1
```

### TM-align 小集合预览

`tmtools==0.3.0` 调用真实 TM-align，适合固定基准或少量本地结构的快速预览：

```powershell
python -m bur01_structmine.cli benchmark --engine tmtools
```

它输出真实 qTM、tTM、RMSD、对齐序列和残基映射。`evalue`、`bits` 和 `prob` 留空，不再伪造 Foldseek 指标。

## 网页工作台

网页包含七个视图：

- **运行**：查询、数据库、阈值、预检、运行、取消、实时事件和错误。
- **NR挖掘**：本机真实命中子集、完整 NR 作业导出/导入、序列证据筛选和送入 AF2 复核。
- **AF预测**：粘贴/上传 FASTA 或选择候选，检查隐私确认与催化位点映射，运行预测或预测并挖掘，查看模型排名、pLDDT、PAE 和三维结构。
- **候选**：按来源、TM-score、覆盖度、WASE 层级、reviewed 状态和证据完整度筛选。
- **结构**：显示查询与候选的 TM-align 叠合，并突出 Ser111/Asp241/His274 及其映射残基。
- **历史**：查看每次运行的状态、命中数、警告、错误和结果目录。
- **诊断**：检查 WSL、Foldseek、AVX2、数据库、网络、磁盘、Python 依赖和写权限。

候选可导出为 Excel、TSV 和 FASTA，也可在界面中重新注释。

### 长期运行方式

- 固定本机地址：`http://127.0.0.1:8501`。
- Windows 登录后自动启动；不要求浏览器或命令窗口始终打开。
- 每 5 秒检查一次健康状态，网页进程退出后自动重新启动。
- Windows 计划任务每分钟触发一次恢复检查；守护程序本身异常退出后最多约 1 分钟恢复。
- 状态和日志保存在 `.cache/ui_service/`，不会进入发布包或覆盖正式结果。
- 默认仅本机可访问。休眠、关机或未登录 Windows 时网页不可访问；恢复并登录后会重新提供服务。

### 在线专家演示

双击 `install_online_ui.bat` 可建立带密码的 HTTPS 演示通道，其他电脑无需安装项目即可实时操作网页。程序会同时打开最新在线网址和本机分享文件 `.cache/remote_ui/SHARE_WITH_EXPERTS.txt`。

默认 Quick Tunnel 的网址在通道重启后可能变化，适合临时展示，不等同于正式云托管。固定网址需要 Cloudflare 命名 Tunnel、自己的域名和 Access 允许名单。操作、安全边界和固定域名步骤见 `ONLINE_DEMO_GUIDE_CN.md`。

## 生产搜索

准备完整 AlphaFold/Swiss-Prot Foldseek 数据库：

```powershell
python prepare_databases.py --config configs/alphafold_swissprot.yaml --only alphafold_swissprot
```

运行生产搜索：

```powershell
python -m bur01_structmine.cli production-run --config configs/alphafold_swissprot.yaml
```

本机内存为 32 GB，因此本项目不下载 UniProt50。默认使用完整 AlphaFold/Swiss-Prot，以避免超出可用内存。

生产配置固定使用官方 AlphaFold/Swiss-Prot `v6 2025_03` 归档，并校验长度和 SHA-256。若 WSL 无法访问 Foldseek 的 Cloudflare 跳转域名，准备器会从配置中的官方 S3 对象断点续传、安全解包并写出来源清单。

## 每次运行的输出

所有正式运行写入：

```text
results/runs/<run_id>/
```

主要文件：

```text
run_manifest.json                 运行清单、输入哈希、软件和数据库版本
config.snapshot.yaml              完整配置快照
events.jsonl                      逐阶段进度事件
*.search.raw.tsv/.xlsx            原始结构命中
*.search.filtered.tsv/.xlsx       阈值过滤结果
*.annotated.tsv/.xlsx             注释与 WASE-screen 结果
*.candidate_proteins.fasta        候选蛋白序列
*.download_ncbi_cds.sh            NCBI CDS 下载辅助脚本
benchmark_report.json             固定基准验收报告，仅基准运行产生
```

`results/runs/latest.json` 指向最近一次完成或失败的运行。失败原因不会被静默吞掉。

## 重要结果字段

结构字段：

- `search_engine`、`metric_source`、`alignment_version`
- `query_chain`、`target_chain`
- `qtmscore`、`ttmscore`、`alntmscore`
- `rmsd`、`qcov`、`tcov`、`alnlen`、`fident`
- `aligned_query`、`alignment_match_line`、`aligned_target`

催化核心字段：

- `mapped_catalytic_residues`
- `mapped_catalytic_identity_count`
- `catalytic_ser_his_distance`
- `catalytic_asp_his_distance`
- `catalytic_geometry_contacts`

注释字段：

- `target_uniprot_accession`、`pdb_id`、`ncbi_protein_accession`
- `uniprot_is_reviewed`、`alphafold_global_plddt`
- `final_protein_sequence`、`final_sequence_source`
- `annotation_status`、`annotation_errors`

## WASE-screen v0.1

WASE-screen 是透明的候选预筛分，不是正式 WASE 模型，也不预测 S/H ratio。

```text
structure_fit_score          45 分
catalytic_core_score         30 分
annotation_confidence_score  15 分
experimental_readiness_score 10 分
```

结构分公式：

```text
25 × min(qTM, tTM)
+ 10 × sqrt(qcov × tcov)
+ 10 × clamp(1 - RMSD / 5, 0, 1)
```

高优先级必须同时满足：

- 总分至少 70。
- 证据覆盖率至少 75%。
- `min(qTM,tTM) >= 0.5`。
- 催化核心至少 20 分。

输出包括：

```text
wase_screen_score
wase_evidence_coverage
wase_screen_uncertainty
wase_tier
wase_missing_evidence
wase_version=0.1.0
```

通道专一性、醇驻留、水排斥、产物逃逸、发酵来源先验和进化证据目前均列为缺失证据，不能据此宣称水相酯化机制已经被验证。

## 注释可靠性

外部注释请求使用：

- 批量 UniProt 请求。
- 对 429 和 5xx 的自动重试。
- 指数退避。
- 30 天磁盘缓存。
- 网络失败时的陈旧缓存回退。
- 明确的 `annotation_status` 和 `annotation_errors`。

网络失败不会再被静默写成空白。

## 验证

```powershell
python -m pytest -q
python -m bur01_structmine.cli doctor --config configs/benchmark_foldseek.yaml
python -m bur01_structmine.cli benchmark --engine foldseek
```

测试使用独立临时目录，不写入正式 `results`。当前验证证据见 `VALIDATION_REPORT.md`。

保留的 v0.3 本机验收（2026-07-11）：

- `43 passed`；`doctor` 为 `31 passed, 0 warnings, 0 failed`。
- 8HFW 标准预测：5 个模型，rank_001 平均 pLDDT `98.07`、pTM `0.950`，相对实验结构 `min(qTM,tTM)=0.9925`。
- Q9FVW3 候选标准预测：5 个模型，平均 pLDDT `83.54`、pTM `0.740`，相对 AlphaFold DB v6 `min(qTM,tTM)=0.7815`，来源运行和 accession 均可追溯。
- 预测 8HFW 固定 Foldseek 基准：最低阳性 `0.8091` 高于最高阴性 `0.2493`，P00520 被拒绝。
- 预测 8HFW 搜索完整 590,183 结构 Swiss-Prot 库：999 个原始命中、12 个默认阈值候选、12 个有效注释，零警告。
- 网页在 `1440x900` 和 `390x844` 下通过截图、横向溢出、键盘焦点、PAE、pLDDT 和 WebGL 非空检查。

v0.4 开发验证（2026-07-24）：

- `65 passed`；无网络 doctor 配置为 `29 passed, 3 expected warnings, 0 failed`。
- MMseqs2 `18.8cc5c` 和 NCBI BLAST+ `2.17.0` 通过真实运行时检查，BLAST+ MD5 与 SHA-256 已记录。
- NR 网页已通过桌面端/移动端截图、横向溢出、键盘焦点、验证数据标识和隐藏临时服务关闭检查；截图位于 `docs/screenshots/nr_workflow_*.png`。
- 真实 NCBI 命中子集等待逐次传输确认，完整 NR 服务器验收仍待执行。
- 开发分支可以推送，但不得合并、创建 `v0.4.0` 标签或 Release。

## 科学边界

- 结构相似不等于相同催化活性。
- WASE-screen 只用于候选优先级，不是活性或选择性预测。
- AlphaFold 结构必须结合 pLDDT、PAE、注释质量和实验验证判断。
- 催化三联体映射是结构证据的一部分，不能替代底物通道、水化和动力学证据。
- 最终候选仍需基因合成、表达纯化、底物谱和水相酯化实验验证。

## 主要文件

- `MODEL_GUIDE_CN.md`：模型用途、原理、创新点和边界的详细中文说明。
- `MODEL_GUIDE_EN.md`：英文模型说明。
- `PROJECT_OVERVIEW.md`：项目结构和入口概览。
- `PROJECT_MANIFEST.md`：完整文件清单。
- `VALIDATION_REPORT.md`：验证结果。
- `COMPLETION_AUDIT.md`：v0.4 完成状态与外部验收门槛。
- `CHANGELOG.md`：版本变更。
- `THIRD_PARTY_LICENSES.md`：第三方软件和许可证。
- `docs/NR_WORKFLOW_CN.md` / `docs/NR_WORKFLOW_EN.md`：NR 操作与证据边界。
- `docs/V0.4_DEVELOPMENT_NOTES_CN.md` / `docs/V0.4_DEVELOPMENT_NOTES_EN.md`：v0.4 当前开发状态。
- `docs/RELEASE_NOTES_CN.md` / `docs/RELEASE_NOTES_EN.md`：v0.3 双语发布说明。
