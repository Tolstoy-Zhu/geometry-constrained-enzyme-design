# Bur01-StructMine V1：PyCharm 可运行代码包

本代码包已经按照 **PyCharm 项目** 重新整理，重点解决 Windows 电脑上打开压缩包后出现乱码的问题。

## 1. 这个代码包能做什么？

本项目用于实现 **Bur01 结构驱动的新酶挖掘 V1 管线**：

```text
Bur01 PDB结构 + 蛋白序列
        ↓
Foldseek 结构相似性搜索 PDB / AlphaFold DB / 自建结构库
        ↓
筛选结构相似蛋白
        ↓
RCSB / UniProt / NCBI 信息反查
        ↓
输出候选蛋白序列、注释信息、NCBI下载命令和候选合成建议
```

注意：

- **PDB数据库** 可以直接做结构搜索；
- **NCBI主要是序列数据库**，不能直接做结构搜索；
- 如果想从 NCBI 全库中做结构相似性挖掘，需要先把候选蛋白序列预测成结构，再用 Foldseek 建库检索；
- V1版本优先实现 **Foldseek结构搜索 + 数据库反查 + 候选输出**。

---

## 2. 文件夹结构

```text
Bur01_StructMine_PyCharm/
├── run_pipeline.py                 # 主入口：完整管线
├── run_demo_no_foldseek.py          # 测试入口：不需要Foldseek和联网
├── config.yaml                      # 参数配置文件
├── requirements.txt                 # Python依赖
├── README_FIRST.txt                 # 快速说明
├── README_CN.md                     # 中文说明
├── data/
│   ├── query/
│   │   ├── Bur01m.pdb               # 示例Bur01结构文件
│   │   └── Bur01m_seed.fasta        # 示例Bur01蛋白序列
│   └── example_foldseek_result.tsv  # 示例Foldseek结果，用于无Foldseek测试
├── src/
│   └── bur01_structmine/
│       ├── pipeline.py              # 主流程
│       ├── foldseek_runner.py       # Foldseek封装
│       ├── results_parser.py        # Foldseek结果解析与过滤
│       ├── rcsb_client.py           # RCSB PDB接口
│       ├── uniprot_client.py        # UniProt接口
│       ├── annotator.py             # 命中结果注释
│       ├── scoring.py               # 候选酶简单评分
│       ├── ncbi_commands.py         # NCBI下载命令生成
│       ├── pdb_utils.py             # PDB序列提取
│       ├── io_utils.py              # 输入输出工具
│       └── config.py                # 配置读取
├── scripts/
│   ├── run_demo_windows.bat         # Windows测试脚本
│   ├── run_demo_linux.sh            # Linux/WSL测试脚本
│   ├── download_pdb100_database.sh  # 下载PDB100数据库
│   └── download_afdb_database.sh    # 下载AlphaFold DB示例脚本
├── databases/                       # Foldseek数据库存放位置
├── results/                         # 输出结果
└── tmp/                             # Foldseek临时文件
```

---

## 3. 在 PyCharm 中打开

### 第一步：解压

将压缩包解压到一个英文路径，例如：

```text
D:\Bur01_StructMine_PyCharm
```

不建议放在含中文或空格很多的路径中。

### 第二步：用 PyCharm 打开

PyCharm 中选择：

```text
File → Open → 选择 Bur01_StructMine_PyCharm 文件夹
```

### 第三步：配置 Python 解释器

建议 Python 版本：

```text
Python 3.10 或 Python 3.11
```

PyCharm 中选择：

```text
Settings → Project → Python Interpreter → Add Interpreter
```

可以选择 Virtualenv。

### 第四步：安装依赖

在 PyCharm 底部 Terminal 中运行：

```bash
pip install -r requirements.txt
```

---

## 4. 不安装 Foldseek 的快速测试

这个测试只检查 Python 环境、PDB读取和文件输出是否正常。

运行：

```bash
python run_demo_no_foldseek.py
```

如果成功，会生成：

```text
results/demo_query_sequence.fasta
```

说明 PyCharm 项目已经可以正常运行。

---

## 5. 用示例 Foldseek 结果测试完整流程

这个测试不需要安装 Foldseek，但会尝试联网访问 RCSB 和 UniProt。

运行：

```bash
python run_pipeline.py --config config.yaml --skip-foldseek --input-tsv data/example_foldseek_result.tsv
```

输出文件：

```text
results/Bur01_vs_database.foldseek.raw.xlsx
results/Bur01_vs_database.foldseek.filtered.xlsx
results/Bur01_vs_database.annotated.xlsx
results/Bur01_vs_database.candidate_proteins.fasta
results/Bur01_vs_database.download_ncbi_cds.sh
```

如果没有联网，RCSB/UniProt注释列可能为空，但前面的解析、过滤和输出流程仍然能跑通。

---

## 6. 正式运行：Linux 或 Windows WSL2 推荐

Foldseek 在 Linux/WSL2 上最稳定。Windows 用户建议安装 WSL2 Ubuntu。

### 6.1 安装 conda 环境

```bash
conda create -n bur01_structmine python=3.10 -y
conda activate bur01_structmine
pip install -r requirements.txt
```

### 6.2 安装 Foldseek

推荐使用 conda 安装：

```bash
conda install -c conda-forge -c bioconda foldseek -y
```

检查是否安装成功：

```bash
foldseek version
```

### 6.3 下载 PDB100 数据库

```bash
bash scripts/download_pdb100_database.sh
```

如果脚本报数据库名不支持，先运行：

```bash
foldseek databases
```

查看当前 Foldseek 版本支持的数据库名称，并修改脚本中的数据库名。

### 6.4 正式搜索

```bash
python run_pipeline.py --config config.yaml
```

---

## 7. 配置文件说明

主要修改 `config.yaml`：

```yaml
query:
  pdb_file: data/query/Bur01m.pdb
  fasta_file: data/query/Bur01m_seed.fasta

foldseek:
  executable: foldseek
  target_db: databases/pdb100
  tmp_dir: tmp
  max_seqs: 1000
  threads: 8

filter:
  min_alntmscore: 0.45
  max_rmsd: 4.0
  min_qcov: 0.60
  min_tcov: 0.40
  min_alnlen: 180
  max_sequence_identity: 0.95
```

如果你要换成自己的 Bur01 文件，把 PDB 放到：

```text
data/query/
```

然后修改：

```yaml
query:
  pdb_file: data/query/your_Bur01.pdb
```

---

## 8. 云算力运行建议

如果只搜索 PDB100，普通电脑或实验室服务器即可。

如果要做 AlphaFold DB 或大规模自建 NCBI 结构库，建议使用云算力。

推荐配置：

| 任务 | 推荐平台 | 配置建议 |
|---|---|---|
| Foldseek vs PDB100 | 本地Linux/WSL2 | CPU 8核，内存16–32 GB |
| Foldseek vs AlphaFold DB | 云服务器 | CPU 16–32核，内存64–128 GB，硬盘500 GB以上 |
| 自建NCBI结构库 | 云服务器 | CPU 32核以上，内存128 GB，硬盘1 TB以上 |
| ColabFold预测结构 | GPU云服务器 | NVIDIA A10/A100/4090，显存24 GB以上 |
| AlphaFold3本地全量部署 | 高性能GPU服务器 | Linux，NVIDIA GPU，数据库空间可能接近TB级 |

V1版本不建议一开始就全量预测 NCBI 蛋白结构，建议先用关键词、HMM或BLAST缩小候选，再预测结构并用 Foldseek 建库。

---

## 9. 输出结果如何使用？

最重要的是：

```text
results/Bur01_vs_database.annotated.xlsx
```

建议重点看这些列：

- `target`：Foldseek命中ID；
- `alntmscore`：结构相似性指标；
- `rmsd`：结构偏差；
- `qcov/tcov`：结构覆盖度；
- `fident`：序列一致性；
- `rcsb_protein_name` / `uniprot_protein_name`：功能注释；
- `organism`：来源物种；
- `final_protein_sequence`：最终蛋白序列；
- `EMBL_accessions` / `RefSeq_accessions` / `NCBI_GeneID`：后续找CDS的入口；
- `synthesis_recommendation`：是否建议进入基因合成。

---

## 10. 获得基因序列的方式

本项目 V1 主要自动整理 NCBI 入口，并生成脚本：

```text
results/Bur01_vs_database.download_ncbi_cds.sh
```

如果候选蛋白有 NCBI GeneID，可用 NCBI Datasets 下载 gene、CDS 和 protein。

如果没有天然 CDS，可根据 `final_protein_sequence` 直接进行密码子优化和人工基因合成。

---

## 11. 常见问题

### 问题1：PyCharm 打开代码乱码

本项目全部文件名使用英文，文本文件使用 UTF-8 编码，并带有 `.editorconfig`。如果仍乱码，在 PyCharm 中设置：

```text
Settings → Editor → File Encodings → Global Encoding = UTF-8
```

### 问题2：运行正式流程提示找不到 Foldseek

说明没有安装 Foldseek，或 Foldseek 不在 PATH 中。

先运行：

```bash
foldseek version
```

如果不行，先安装 Foldseek。

### 问题3：RCSB/UniProt 注释为空

可能是没有联网，或目标ID不是PDB/UniProt格式。先检查网络，再检查 `target` 列格式。

### 问题4：NCBI CDS 不能直接获得

蛋白序列不一定都有清晰的天然基因记录。对于实验表达，最稳妥的方法是：

```text
完整蛋白序列 → 表达宿主密码子优化 → 人工合成基因
```

---

## 12. V1版本边界

已实现：

- Bur01 PDB输入；
- Foldseek结构搜索封装；
- Foldseek结果解析和过滤；
- RCSB/UniProt注释反查；
- NCBI GeneID/RefSeq入口整理；
- 候选蛋白FASTA输出；
- 候选合成推荐；
- NCBI下载命令生成。

尚未完全自动化：

- NCBI全库蛋白结构预测；
- AlphaFold/ColabFold批量预测结构；
- 真实CDS的一键下载和清洗；
- 催化三联体的三维空间自动识别；
- 与DALI结果的一键融合。

这些可以作为 V2 版本继续开发。
