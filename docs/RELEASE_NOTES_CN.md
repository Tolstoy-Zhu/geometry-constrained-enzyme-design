# AquaEsterNet-StructMine v0.3.0 发布说明

## 本版定位

v0.3.0 在 v0.2 可追溯结构挖掘流程上增加本地 AlphaFold2-PTM 单体预测。用户可以输入一条 FASTA 或选择已有候选，执行“仅预测”或“预测并挖掘”，并把 rank_001 结构作为同一 Foldseek 服务层的动态查询。

## 主要变化

- LocalColabFold 1.6.1、AlphaFold2-PTM、JAX/CUDA 和 OpenMM 在 `Ubuntu-Foldseek` 中隔离运行。
- 输入限制为单条 20-1000 aa 蛋白序列，`X` 不超过 5%；格式错误一次性列出。
- 标准档固定为 5 个模型、3 次 recycle、1 个 seed、rank_001 Amber 松弛；模板默认关闭。
- 公共 MSA 与本地 GPU 推理解耦；MSA 按序列哈希缓存、严格串行、逐任务确认第三方传输，并最多三次指数退避。
- 每台 GPU 只运行一个预测任务；支持取消、失效锁恢复和明确的显存不足错误，不回退 CPU。
- 自动将 8HFW `Ser111/Asp241/His274` 映射到新序列，允许人工修正；身份不符时不伪造催化证据。
- 预测与 Foldseek 运行清单双向关联，记录软件版本、权重来源、SHA-256、模型排名、pLDDT、pTM、PAE、事件、警告和失败原因。
- 网页新增“AF预测”视图，以及候选重新预测、模型排名、pLDDT 曲线、PAE 热图和按 pLDDT 着色的三维结构。
- 新增 Python 3.11 非 GPU CI、假 ColabFold 端到端测试、脱敏发布包构建和 Playwright UI 验收。

## 本机验证

- 自动测试：34 项通过。
- 环境诊断：31 项通过，0 警告，0 失败。
- 8HFW：5 模型，rank_001 平均 pLDDT 98.07，pTM 0.950，相对实验结构最低双向 TM-score 0.9925。
- Q9FVW3：5 模型，rank_001 平均 pLDDT 83.54，pTM 0.740，相对 AlphaFold DB v6 最低双向 TM-score 0.7815。
- 预测 8HFW 基准：最低阳性 0.8091，高于最高阴性 0.2493；P00520 被拒绝。
- 完整 AlphaFold/Swiss-Prot：590,183 个结构，999 个原始命中，12 个过滤候选，12 个有效注释，0 警告。
- 网页：1440x900 与 390x844 均通过截图、横向溢出、键盘焦点、图表和 WebGL 非空检查。

## 隐私与许可

缓存未命中的序列只有在用户逐次确认后才发送到公共 ColabFold MSA 服务。AlphaFold 参数使用随附的 CC BY-NC-SA 4.0 许可，不包含在仓库或 Release 附件中。发布包同样排除数据库、缓存、Obsidian 文件、真实候选结果、模型权重、运行日志和本机路径。

## 科学边界

AlphaFold 预测和结构相似性不证明酶活性、S/H ratio 或水相酯化机制。WASE-screen v0.1 仍是透明预筛分，不是正式 WASE 模型。候选必须结合 pLDDT、PAE、结构覆盖、注释质量和实验验证解释。

## 升级说明

v0.2 结果、数据库、Obsidian 库和既有交付包保持不变。v0.3 的预测写入 `results/predictions/`，后续结构挖掘仍写入 `results/runs/`。
