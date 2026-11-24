# WASE-screen 0.2.0 评分说明

适用程序：AquaEsterNet-StructMine 0.5.0.dev0。分数用于候选预筛排序，不预测酶活、S/H 或实验成功概率。

## 计算依据

| 项目 | 分值 | 依据 |
|---|---:|---|
| 结构匹配 | 45 | 25×min(qTM,tTM) + 10×sqrt(qcov×tcov) + 10×clamp(1−RMSD/5,0,1) |
| 催化核心 | 30 | 映射 Ser 中心 GXSXG 模体 5；映射 Ser/Asp/His 身份各 5；两组几何距离各 5 |
| 注释可靠性 | 15 | 相关功能 5；reviewed/实验结构 5；结构可信度 5 |
| 实验准备度 | 10 | 序列可用且无 X 4；长度 180–500 aa 2；tcov≥0.8 2；可追溯 ID 2 |

输入 TM 和覆盖度进入计算时限制在 0–1；原始搜索字段保留，不改写。Foldseek 原始 `alntmscore` 可能因近似算法略大于 1，展示双向 qTM/tTM 时保留各自定义。

高优先级同时要求总分≥70、评分字段可评估权重≥75、min(qTM,tTM)≥0.5、催化核心≥20。得分≥50但未满足高优先级条件为中优先级；可评估权重<60时首先标为证据不完整。结构筛选仍单独遵循运行配置中的 TM、RMSD、覆盖度、对齐长度和一致性阈值。高优先级不等于通过所有结构筛选。

## 模体和催化几何

TM-align 将查询位点映射到目标链，模体必须来自该 Ser 两侧各两个连续的结构残基；未解析缺口、窗口不足或无法映射记为缺失。全序列 `has_GXSXG_like_motif` 仍保留作辅助信息，不直接得催化分。

Ser 的 OG、Asp 的 OD1/OD2、His 的 ND1/NE2 均需按对应距离所需原子检查；残基身份不符不冒充原子缺失。Ser-His 和 Asp-His 各使用相关原子的最短距离，2–5 Å 为当前筛选范围。此范围只表示几何接近，不证明质子传递或氢键。

- `compatible`：两项均可计算并符合范围。
- `out_of_range`：两项均可计算，至少一项不符合。
- `not_assessable`：无法完整计算；总接触数字段为空。若一项可计算，仅该项计分和计入可评估权重，缺失的另一项列出原因。

Foldseek 原始指标不被 TM-align 重写；补充位点映射记录 `catalytic_mapping_source` 和结构 SHA-256。只有本地坐标可用时才补算；重新评分不自动下载结构或启动 Foldseek。

## 缺失程度与预测误差

`wase_information_missingness = 100 − wase_evidence_coverage`。这里的 coverage 是上述100分项目中可评估字段的权重，不是全部生物学证据的完整程度。通道、水化、乙醇驻留、底物相容性、生态和进化等尚未评估的证据仍列在 `wase_missing_evidence` 中，即使 coverage 为100也不消失。

新记录的 `wase_screen_uncertainty` 和 `wase_prediction_error` 均为空，`wase_prediction_error_status=not_estimated`。旧字段仅在界面兼容显示为“信息缺失程度”。没有统计校准，不提供成功概率或预测误差区间。

## 历史结果

候选页点击“另存新版评分”，或运行：

```powershell
python -m bur01_structmine.cli rescore --run-dir results/runs/<run_id>
```

生成 `results/rescores/<独立编号>/rescored.tsv`、Excel 和 `rescore_manifest.json`。包含原运行编号、原评分版本/分数、分数变化和输入哈希。旧表、旧清单不修改；重新注释也生成新的运行。旧预测缺少催化参考映射时，新记录明确标记缺失，不默认套用8HFW编号。
