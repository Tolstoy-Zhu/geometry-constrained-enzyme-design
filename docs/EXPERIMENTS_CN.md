# 实验模板与评估

网页进入“实验评估”下载 Excel 模板，也可使用 `data/templates/experiment_template.xlsx`。模板首个工作表 `measurements` 记录数据，`units` 列出允许单位。支持 `.xlsx` 和 UTF-8 CSV；单文件上限20 MiB。不支持含宏的工作簿。

## 每行的含义

每行是一条测量记录。列名采用稳定的英文标识，界面默认中文。

| 字段 | 含义 |
|---|---|
| sample_id | 样品标识，必填；保留前导零 |
| candidate_id | 关联候选 accession/ID；未关联时留空并列出提示 |
| metric / value / unit | 指标、有限数值、单位，必填 |
| conditions | 条件，必填；写全底物、温度、pH等；可用JSON对象 |
| method | 测量方法，必填；不同方法不合并 |
| replicate_id | 重复编号，必填；同组重复编号不能重复 |
| predicted_value / prediction_source | 可选预测值及来源；有数值但无来源时不配对 |
| predicted_metric / predicted_unit | 可选，默认与实测相同；不一致时不配对 |
| wase_screen_score | 可选0–100预筛分，仅用于排序比较 |
| wase_source_run_id / wase_version | 可选分数来源运行编号与评分版本；网页关联时记录 |
| data_kind | `experimental` 或 `simulated`，默认experimental；示例强制为simulated |

不要把 WASE-screen 填作活性预测。网页可按精确的候选ID补充当前候选表的 WASE-screen 分数，不覆盖已填数值；含歧义的候选不自动匹配。

## 指标和单位

| metric | 可选单位 |
|---|---|
| activity | U/mL、U/mg、U、kat、kat/kg、umol/min、mmol/L/h |
| relative_activity | %、fraction |
| Tm / optimal_temperature | C、K |
| residual_activity | %、fraction |
| solubility | %、fraction、mg/mL、g/L |
| S/H | ratio |

中文“活性、相对活性、最适温度、残余活性、可溶性”和Topt可作为指标别名。℃/°C规范为C；不自动做单位换算。Tm 与最适温度是不同指标。`fraction` 的1表示100%。未知条件不能可靠合并，需填写为明确且独立的条件标识。

## 汇总和比较

- 同一样品、候选、指标、单位、条件、方法和数据类别的重复，输出均值、样本标准差（ddof=1）、n与原始行号。单次测量标准差为空；重复测量不当作多个独立样品。
- 条件字符串除多余空格外保持原样；JSON对象按键排序。不会猜测“pH7”和“pH=7”是否相同。
- 预测与实测按同一指标、单位、条件、方法和预测来源比较。重复中的冲突预测不平均处理，列出原因。
- MAE=mean(abs(y−p))；RMSE=sqrt(mean((y−p)²))；R²=1−SSE/SST，不是Pearson²，可能为负。
- Pearson和Spearman（并列使用平均秩）只在至少3个有效独立样品且两侧数值均不恒定时给出。此限制同样用于R²；小样本仍可报告MAE/RMSE，但不能据此宣称泛化能力。
- WASE-screen 与活性类指标或S/H仅报告Spearman及Top-k，不报告MAE/RMSE/R²，也不推断成功概率。
- 先选择一个指标/单位/条件/方法组，再设置成功阈值和Top-k。预测排序按成功方向，WASE仍按分数降序；k超过样本数时使用有效样本数。第k名存在跨分界并列时，命中率留空并提示。

## 错误和导出

表格显示工作表/CSV记录对应的行号（含表头）和错误列；错误行不参加统计，但保留在raw页。重复编号的两条记录都排除。无法配对原因见unmatched页；常见原因包括无预测、无预测来源、指标/单位不一致、重复预测冲突、未关联候选。

完整Excel包含raw、valid、errors、aggregates、comparisons、unmatched和summary；另可导出选中分组与阈值。导出文本采取Excel公式转义，不执行上传内容。当前输出不是自动写入生产候选库。

```powershell
python -m bur01_structmine.cli evaluate-experiments --input data/examples/simulated/experimental_measurements.csv --output results/evaluations/example_01 --top-k 3 --success-threshold 30
```

目录已存在时拒绝覆盖。导入错误时仍输出审查文件，并返回退出码2。跨多个条件组时不能共用一个阈值，需先拆分输入。

`data/examples/simulated/` 仅含模拟示例，与真实实验分开；未从论文或历史研究中伪造本平台的验证数据。
