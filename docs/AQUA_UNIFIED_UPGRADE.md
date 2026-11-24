# AquaEsterNet 0.6 工程预览版

本分支在原 StructureMine 工作台中增加项目、家族、远端任务、模型评估和机制计算。软件功能、预测准确性、机制采样分别验收。旧报告与原始结果目录保留，不用新版结果覆盖旧结果。

## 已接入的功能

- SQLite 项目、酶及表达构建、反应、实测记录、任务与模型索引。记录冲突会拒绝覆盖；检测限与未测样品单独保留。
- Bur01、CALB（1TCA A）及 GH1 β-葡萄糖苷酶（1QOX A）配置。GH1 使用 E166/E355；WASE-screen 仅用于 Bur01。
- 网页、命令行共用持久任务服务。SSH 上传明确选定的输入，关闭浏览器不终止计算；重连先核对任务编号，机制计算与特征任务保留同一任务的进度。
- SWI、TmProt、EnzRank 保持独立。ProteinMPNN/LigandMPNN 增加完整组合序列评分；此评分描述序列相容性，不代表已测得的突变上位性。
- CataPro 数据冻结、单位和化合物标准化、MMseqs2 同源分组、冻结 ESM-2 特征、独立 kcat/Km 的 Ridge 与 ExtraTrees 训练、校准、分组评估及模型登记。
- SABIO-RK 新版 export API 的分页 JSON/SBML 下载，以及需要确认构建的参数映射。旧 sabioRestWebServices 地址已停止返回预期数据。
- 已准备 Amber 系统的 OpenMM 8.2.0 分段生产模拟、环境匹配的检查点恢复、逐重复水/醇占据、驻留及几何通道分析。

## 运行与服务器部署

本机沿用原网页命令：

```powershell
python -m bur01_structmine.cli web --port 8501
python -m bur01_structmine.cli families
python -m bur01_structmine.cli projects index --root "<原仓库目录>"
```

CPU 训练依赖单独安装 `requirements-learning.txt`。通过 `AQUA_LEARNING_PYTHON` 指定该环境；`AQUA_PROPERTY_PYTHON` 指向既有性质模型环境。ESM 特征可通过 `AQUA_FEATURE_PYTHON` 指定 GPU 环境。

Windows 多环境部署可复制 `configs/unified/runtime.example.json` 为 `.cache/unified/runtime.json`，填写本机解释器和模型位置后双击 `start_unified.bat`。它启动同一工作台，只监听本机地址；不会改动系统环境变量或替换已有开机任务。端口属于其他仓库时拒绝覆盖，可用 `scripts/unified/start_workbench.ps1 -Port 8502` 指定另一端口。

`data/families/examples` 提供三个公开结构链的序列、项目导入 JSON 和结构检查示例。可以运行 `python scripts/unified/public_examples.py --output <目录>` 重新生成。这些是结构中解析出的序列，不代替完整表达构建记录，不包含虚构的实验测量。

远端使用数据盘中的独立部署目录：

1. 将公开仓库放在 `<部署目录>/code`。
2. 执行 `bash code/scripts/unified/bootstrap_autodl.sh <部署目录>`。
3. 将 `scripts/unified/aqua-python` 复制为 `<部署目录>/bin/aqua-python` 并赋予执行权限。该入口设置模型目录、工作目录及 CUDA 库路径。
4. 模型文件按 `property_assets.py` 和 `design_models.py` 中的来源与校验值准备，分别放在 `models/properties-v1` 与 `models/design-v1`。
5. 在网页“项目与任务”保存 SSH 设置，远端 Python 字段使用上述 `bin/aqua-python`。配置写入本机 `.cache/unified/servers.json`，不提交 Git。
6. 运行 `bin/aqua-python code/scripts/unified/gpu_smoke.py --output logs/gpu-smoke.json`。这是数值环境和检查点检查，不是酶轨迹验收。

可用 `PIP_INDEX_URL`、`AQUA_TORCH_INDEX`、`PIP_FIND_LINKS` 指定可信镜像或已校验的本地安装包。不会修改 GPU 驱动。所有重计算任务默认共用单个 GPU 执行锁。

新项目配置服务器后默认使用远端；旧命令默认本地。run/predict/design/properties 接受 `--family`、`--executor`、`--reaction-config`。带远端或反应配置时生成持久任务。大型已有远端输入用 `remote:/部署目录内的绝对路径`，不会再次上传。任务目录外的远端路径会拒绝。

## 数据与模型

CataPro 原始包：
https://zenodo.org/records/14894710

- 发布 MD5：`2f8411d35fba9f5a9b3069dac59d3010`
- 本次下载 SHA-256：`aa58dd9730fff51e510ea0738af00cbf256aa7cd55d9cecfdcb95932d23ed519`
- 原表没有完整测定条件；只用于无条件基线。

```text
python -m bur01_structmine.cli dataset fetch-catapro --output <原始zip>
python -m bur01_structmine.cli dataset build --spec <数据集JSON> --executor autodl
python -m bur01_structmine.cli dataset features --spec <特征JSON> --executor autodl
python -m bur01_structmine.cli train submit --spec <训练JSON> --executor autodl
python -m bur01_structmine.cli tasks sync --task <任务编号>
```

数据集 JSON 输入为 `{"input":"原始zip路径","source":"catapro","endpoint":"kcat","mmseqs":"mmseqs"}`；Km 独立构建。特征输入为 `{"dataset":"数据集目录","device":"cuda"}`；训练输入为 `{"dataset":"数据集目录","features":"特征目录","conditions":false}`。

MMseqs2 18-8cc5c 冻结为 40% 一致性、双向 80% 覆盖度，全对全搜索边形成连通分组；同一亲本与重复测量继续并组。分组比例 70/10/10/10，数量取整；记录实际组数和样本比例。使用训练/验证选择参数，校准集只计算区间，测试集评估一次。在拟合前按固定随机种子保留 10% 底物身份，这些底物的测量不用于拟合、参数选择或区间校准。训练组中已见序列与保留底物的组合构成独立“新底物”测试，主同源测试中的陌生底物构成“酶和底物同时陌生”测试。主同源分组保持不变；报告实际参与各部分的样本数，并保存底物保留清单。

独立端点发布；只有 ExtraTrees 在验证集入选、封存测试的总体 MAE 低于 Ridge，且按序列/亲本组重采样的 MAE 改善区间支持优势时，登记才会替换默认版本。其他模型标为 experimental。酯合成初速、水解初速及合成/水解比无足够标签，保持尚未训练。

SABIO-RK 使用 [当前 OpenAPI](https://sabiork.h-its.org/openapi/export-api.json)：

```text
python -m bur01_structmine.cli dataset fetch-sabio --query 'EntryID:123' --output <快照目录>
python -m bur01_structmine.cli dataset map-sabio --input <快照目录> --spec <构建映射JSON> --output <映射目录>
```

映射以 SABIO 条目 ID 为键，必须给出 sequence、sequence_source、construct_confirmed；突变体还需要 mutations_confirmed。多底物条目无法唯一归属时，使用 parameter_substrates 将参数名指向 SABIO 化合物 ID。范围值和未映射项留在待核对表。合并不同来源前，须用原论文、构建、底物、条件和测量编号核查重复；目前不自动把来源不明的两行合并。

## 机制计算约束

预算已经按用户决定改为 **120 ns**：亲本 +3 条候选，每个 3 个独立重复，各 10 ns；平衡另计。运行前必须冻结实际完整构建序列、三条候选选择依据、反应条件及独立随机种子。

公共接口要求 prepared.json 包含实际 prmtop、平衡后 inpcrd、同一坐标的 protein_pdb、反应及各文件 SHA-256。构建、装配、参数、几何、平衡五项证据必须通过，未解决项必须为空。会核对 Amber 体系的实际蛋白序列及重原子坐标；旧记录中的 PASS 字样不直接放行。

准备适配器由部署配置 `AQUA_PREPARATION_ADAPTER` 指定，含脚本路径、SHA-256 和 Python。私有研究脚本与原始数据不随公共仓库发布。未验证的金属、共价配体、非标准残基会停止准备。

按不超过 10 ns 分段保存检查点。硬件、OpenMM 或运行环境变化时拒绝直接加载旧二进制检查点；可移植状态继续运行需要单独记录新的运行身份。分析前处理周期边界；多链体系每个活性中心分别选择区域。输出模拟描述量，未进行化学反应能垒预测。

几何通道算法为网格最宽路径，需另行检查网格分辨率、探针半径敏感性。MSM/TPT 必须独立检查重复数、转移事件、隐含时间尺度和 Chapman–Kolmogorov；当前预览不会生成未经检验的速率。

## 验收记录

- 原固定性质模型环境重算 24 条设计及亲本：25 条蛋白、50 组酶—底物；SWI、TmProt、EnzRank 与原报告数值一致。
- 单元测试覆盖编号、构建映射、家族规则、单位、检测限、模拟标签排除、并发任务、重试、取消、结果传输与原文件保留。
- 2026-11-24 已完成 GPU 性质预测、ESM 特征中断恢复、两个自主模型训练和 CALB/GH1 组合设计示例。机器地址、私有构建和候选结果留在研究工作区；公开汇总见 [本次验收记录](validation/2026-11-24-engineering-preview.json)。
- 172 项非 GPU 测试及独立数值环境中的 8 项测试通过。主测试环境未安装 RDKit，其学习模块由独立数值环境运行。
- 原 24 条核心设计映射至 294 残基完整构建后，完成亲本加 24 条候选的性质及 kcat/Km 预测。新结果单独保存；其序列适用域尚未建立，不发布校准区间。
- OpenMM CUDA 与二进制检查点小体系测试通过。120 ns 酶轨迹尚未开始：完整构建坐标、独立结构复核、配体/缓冲液参数及催化几何和平衡检查仍需完成。
- 完整 NR 数据库验收独立保留；当前小规模示例通过不等于完整 NR 验收。

### 本次冻结模型评估

数值均为 log10 尺度；该比较针对冻结数据上的 Ridge 对照，不代表已验证 Bur01 的实验活性。

| 指标 | 封存测试样本 | Ridge MAE | ExtraTrees MAE | ExtraTrees RMSE | Spearman | 按组 MAE 差值 95% 区间 |
|---|---:|---:|---:|---:|---:|---|
| kcat | 2,369 | 1.620 | 1.114 | 1.446 | 0.304 | −0.484 至 −0.226 |
| Km | 3,887 | 1.011 | 0.812 | 1.044 | 0.532 | −0.306 至 −0.153 |

新底物、酶与底物同时陌生的结果、数据划分数量及模型版本分别保存在公开验收 JSON 中。原始公共数据仅具整理数据行级来源；原始论文测定条件和跨库重复测量仍需逐项映射。条件增强模型和机制特征增强模型尚无对应训练验收。
