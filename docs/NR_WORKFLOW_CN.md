# AquaEsterNet-StructMine v0.4 NR 使用说明

## 1. 这项功能解决什么问题

结构库能直接发现与 Bur01 折叠相近的蛋白，但许多 NCBI NR 条目只有序列，没有实验结构或 AlphaFold 结构。v0.4 增加独立的 NR 序列扩增层，用已确认的结构候选作为种子，在更大的天然序列空间中寻找远缘候选。

NR 层只回答三个问题：

1. 哪些天然序列与 Bur01 或结构种子具有可追溯的序列匹配？
2. 哪些代表序列保留 GXSXG-like 模体及映射后的 Ser111/Asp241/His274？
3. 哪些序列值得进入 AF2 和结构复核？

它不直接回答“结构是否相似”“是否具有水相酯合成活性”或“S/H ratio 是否提高”。这些结论仍需要 AF2、Foldseek/TM-align、WASE-screen 和实验验证。

## 2. 五个证据层必须分开

| 层 | 输入与输出 | 可以说明 | 不能说明 |
|---|---|---|---|
| 结构挖掘 | PDB/AlphaFold + Foldseek/TM-align | 折叠、覆盖度、RMSD、催化核心映射 | 活性和反应方向 |
| NR 序列扩增 | 结构种子 + MMseqs2/NCBI BLAST | 序列邻近、聚类、模体、三联体身份 | 结构相似或 WASE 高优先级 |
| NR 新颖性审计 | 预测/设计序列 + NR 最近邻 | identity、coverage、E-value、最近邻来源 | 天然序列自身命中的创新性 |
| AF2 预测 | 单体 FASTA + 公共 MSA + 本地 GPU | 结构假设、pLDDT、pTM、PAE | 实验结构或催化活性 |
| WASE-screen | 已复核结构和注释 | 透明候选优先级与缺失证据 | 正式 WASE 或 S/H ratio 预测 |

所有未经结构复核的 NR 候选均输出：

```text
wase_tier = Not evaluated - structure required
wase_scope = not_evaluated_sequence_only
```

## 3. 复现预设

v0.4 的 `hong2025` 预设参考 VenusMine 的公开 MMseqs2 搜索与聚类步骤：

```text
MMseqs2 iterations       2
max-seqs                 2500
search-type              1
cluster min sequence id  0.50
```

当前实现不引入未经训练验证的 PLM 综合分数。序列排序只使用 E-value、bits、identity、query/target coverage、长度、GXSXG-like 模体、催化三联体身份和注释状态。

参考资料：

- VenusMine Step 3: https://github.com/ai4protein/VenusMine/tree/main/Step3%20MMseqs%20Search
- VenusMine Step 4: https://github.com/ai4protein/VenusMine/tree/main/Step4%20MMseqs%20Cluster
- VenusMine 论文: https://doi.org/10.1038/s41467-025-61599-z
- CPDiffusion: https://doi.org/10.1038/s41421-024-00728-2
- MMseqs2: https://github.com/soedinglab/MMseqs2
- NCBI BLAST API 规范: https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html

## 4. 本机演示模式

本机不下载完整 NR。`nr-demo-build` 对 Bur01 执行一次真实 NCBI 远程 `blastp nr`，最多保留 2,500 个命中，批量取回序列，并写成可重复使用的本地快照。

远程请求必须同时满足：

1. 提供有效的 `AQUA_NCBI_EMAIL` 或 `--email`。
2. 每次缓存未命中时显式使用 `--confirm-ncbi`。
3. 请求严格串行，轮询间隔至少 60 秒，失败采用有限重试。

```powershell
$env:AQUA_NCBI_EMAIL = "name@example.org"
python -m bur01_structmine.cli nr-demo-build `
  --query-fasta data/query/Bur01m_seed.fasta `
  --max-hits 2500 `
  --confirm-ncbi
```

快照写入：

```text
databases/nr_demo/<snapshot_id>/
  nr_remote_snapshot.fasta
  nr_remote_hits.tsv
  nr_snapshot_manifest.json
databases/nr_demo/latest.json
```

清单记录查询哈希、RID、参数、创建时间、序列数、文件 SHA-256 和明确的 `nr_remote_snapshot` 范围。它不是完整 NR。

从已完成的结构运行扩增：

```powershell
python -m bur01_structmine.cli nr-expand `
  --source-run <structure_run_id> `
  --config configs/nr_local_demo.yaml
```

## 5. 完整 NR 服务器模式

完整 NR 只通过显式作业包运行。本机不会在普通搜索时自动下载 NR，也不保存 SSH 密码、令牌、主机地址或私钥。

导出：

```powershell
python -m bur01_structmine.cli nr-export-job `
  --source-run <structure_run_id> `
  --config configs/nr_full_server.yaml
```

作业包包含：

```text
seeds.fasta
bur01_query.fasta
config.snapshot.yaml
job_metadata.json
resource_precheck.sh
prepare_nr_database.sh
run_nr_pipeline.sh
submit.slurm
submit.pbs
nr_server_postprocess.py
BUNDLE_SHA256SUMS
README_CN_EN.md
```

服务器首次建库时显式运行：

```bash
bash prepare_nr_database.sh
```

该脚本调用：

```bash
mmseqs databases NR "$NR_DB" "$TMP_DIR"
```

然后运行通用 Bash、Slurm 或 PBS 入口。`NR_SNAPSHOT_DATE` 是强制字段。服务器结果包含 `server_result_manifest.json` 和 `RESULT_SHA256SUMS`。

导回：

```powershell
python -m bur01_structmine.cli nr-import-job `
  --bundle <nr_run_id>_server_result.zip
```

导入器拒绝路径穿越、缺失文件、SHA-256 不一致、非 `full_nr` 范围、错误版本和错误 run ID。正式 `v0.4.0` Release 只有在完整 NR 结果导回并通过验收后才能创建。

## 6. 序列处理原理

1. 从结构运行导出 Bur01 查询序列及所有可用 Foldseek 目标序列。
2. 先按 accession 选择最佳命中，再按完整序列 SHA-256 去重。
3. 以 50% 序列一致性聚类，保留代表序列和完整成员映射。
4. 每个代表序列重新与 8HFW chain A 直接对齐。
5. 分别映射 Ser111、Asp241、His274 的目标位置、观察残基和缺失原因。
6. 按透明阈值生成 `AF2-ready` 队列。

P00520 等远缘负对照必须因 E-value、覆盖度或催化三联体证据不足而停留在 `hold`，不得通过默认序列预筛。

## 7. AF2 复核

网页“NR挖掘”视图可把代表序列送入“AF预测”。命令行：

```powershell
python -m bur01_structmine.cli predict `
  --nr-run <nr_run_id> `
  --accession <accession> `
  --preset standard `
  --confirm-public-msa
```

预测清单记录 `source_nr_run_id` 和 accession，NR 清单反向记录 `child_prediction_ids`。使用 `predict-and-mine` 时，rank_001 再进入 Foldseek/TM-align 和 WASE-screen，形成 NR → AF2 → 结构运行的双向谱系。

## 8. 新颖性审计

```powershell
python -m bur01_structmine.cli nr-novelty `
  --fasta designed_sequence.fasta `
  --backend ncbi_remote `
  --email name@example.org `
  --confirm-ncbi
```

输出最近邻 accession、identity、coverage、E-value、标题/物种和数据库快照。天然 NR 来源序列出现自身 100% 命中是正常情况，不得作为“创新性”证据。

## 9. 输出目录

每次扩增写入 `results/nr_runs/<nr_run_id>/`：

```text
nr_run_manifest.json
config.snapshot.yaml
seeds.fasta / seeds.tsv
nr_hits.raw.tsv
nr_hits.deduplicated.tsv
nr_s50_rep_seq.fasta
nr_s50_cluster.tsv
bur01_vs_representatives.tsv
nr_candidates.sequence_screen.tsv/.xlsx
nr_candidates.fasta
af2_ready.tsv
events.jsonl
SHA256SUMS
```

数据库、真实候选、缓存、日志、模型权重、服务器配置和 ZIP 均由 `.gitignore` 排除，不进入 GitHub 或公开 Release。

## 10. 当前验收状态

截至 2026-07-24：

- MMseqs2 `18.8cc5c`：已安装并通过诊断。
- NCBI BLAST+ `2.17.0`：已按官方 MD5 校验并记录 SHA-256。
- 自动测试：65 项通过；覆盖命令、缓存、限速、去重、聚类、催化映射、清单、作业包校验、新颖性审计和 NR→AF2 谱系。
- 本机真实 NCBI 快照：需要用户逐次确认并提供联系邮箱后运行。
- 完整 NR：等待服务器/HPC 结果导回；因此当前只允许开发分支，不合并、不打 `v0.4.0` 标签。

`configs/ncbi_prostt5.yaml` 仅保留为小型自定义 FASTA 的 legacy 入口，不代表完整 NR。
