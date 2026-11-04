# AquaEsterNet 本地酶序列设计

AquaDesign 针对具有 Bur01 折叠和 Ser–Asp–His 催化核心的酶，提出适合逐个实验比较的少量替换。AquaEsterNet 管理位点选择、同源约束、替换比较和坐标复核，使用公开的 ProteinMPNN 与 LigandMPNN 模型做局部序列推理。设计时序列和结构均在本机处理。[模型代码与许可](https://github.com/dauparas/LigandMPNN/tree/26ec57ac976ade5379920dbd43c7f97a91cf82de) · [Bur01 参考结构](https://www.rcsb.org/structure/8HFW)

## 安装和打开

在项目的 Python 环境安装设计依赖，然后获取模型。命令只在安装阶段联网。两种模型权重及上游推理代码按 SHA-256 校验，默认保存在项目的 .tools/aqua-design，可用 AQUA_DESIGN_MODELS 指定其他目录。本机已验证 Python 3.11、PyTorch 2.9.1 CPU。已有 CUDA 版 PyTorch 的环境可以通过命令行显式选择 CUDA。

~~~sh
python -m pip install -e ".[design]"
python -m bur01_structmine.cli design-setup
~~~

打开完整工作台，在“结构”页展开“本地酶序列设计”，即可选用公开 Bur01、候选表中有本地坐标和完整序列的亲本，或上传 PDB/mmCIF。选择蛋白链，可再选择带有实际原子坐标的配体。运行后可下载设计 FASTA、位点表和全部替换分数。在“AF预测”页将输入方式切换为“本地设计序列”，选择批次和候选，沿用原预测流程与公共 MSA 的确认机制。

## 设计方法

1. 从第一模型中提取所选单链的实际坐标序列。与 8HFW 链 A 成对 TM-align，要求双向 TM-score 的较小值至少为 0.5，并逐一定位 Ser111、Asp241 和 His274。缺少必要主链原子、催化侧链原子或出现链断裂时停止。独立上传的坐标按其已解析链段作为设计构建，需与拟表达的构建核对。来自候选表的结构另与候选序列作精确比对。
2. 固定催化三联体、映射到的 Ser 中心 GXSXG 模体，以及任一重原子距催化侧链功能原子不超过 4 Å 的残基。可指定额外固定位置。无配体时按到催化核心的距离选择位点，有配体时按到所选配体重原子的距离选择，默认在 10 Å 内最多评价 16 个位点。距离取自亲本的一份静态结构。
3. 可提供以亲本为第一条的对齐 FASTA 或 A3M。移除小写插入，按亲本有效列裁齐，去除完全重复序列和亲本。某列至少有 3 条非亲本同源序列提供标准氨基酸、覆盖比例至少 70%、且至少 85% 的有效残基与亲本相同，才将该列列入保守保护。没有比对时只用结构和人工约束。上述阈值是筛选设置。
4. 在其余亲本残基已知的条件下，把待评位置放在最后解码。分别计算各替换相对亲本残基的对数概率差。无配体时用 ProteinMPNN。有配体时用 LigandMPNN 在相同权重、骨架和解码顺序下分别保留、遮蔽配体原子，两次对数概率差的差值描述模型对当前配体坐标的敏感程度。单独的化学名称或 SMILES 没有三维坐标，无法提供此项输入。
5. 默认每个位点保留相容性最高的非亲本残基。选“配体遮蔽对照”时，先取有配体条件下相容性最好的三种替换，再选其中配体差值最高者。可生成 1–3 处替换的设计。组合只使用亲本 Cα 间距至少 8 Å 的位点，按各自突变数分组排序、轮流取样。多点分数采用单点贡献之和，未估计上位性。候选不足时输出实际数量。

## 示例与产物

~~~sh
python -m bur01_structmine.cli design --parent data/query/Bur01m.pdb --chain A --fixed A240 --max-edits 2 --count 12 --output results/designs/bur01_example
~~~

有带配体的复合物时，可追加 --ligand B:ETH:1 或多次指定 --ligand，并使用 --priority ligand_response。配体标识必须与文件的链、残基名、编号及插入码一致。可用 --msa homologs.a3m 引入家族信息，用 --seed 固定解码顺序。CLI 默认使用 CPU。

| 文件 | 可回答的问题 |
|---|---|
| sites.tsv | 每个位置为何固定、进入模型或排在候选位点上限之外？ |
| substitutions.tsv | 各位点的 19 种替换各有多大的模型对数概率差？显式遮蔽配体后相差多少？ |
| candidates.tsv、candidates.fasta | 哪些序列入选？具体替换对应亲本结构和序列的哪些编号？ |
| manifest.json | 使用了哪条链、哪些配体、哪些固定位置、哪个权重和随机种子？ |
| reviews/ | 哪一份坐标与该设计序列完全匹配？和亲本的双向 TM-score、RMSD、Ser–His 与 Asp–His 侧链距离是多少？ |

候选表同时给出 GRAVY 和估计的 pH 7 净电荷，用于描述序列性质。配体差值受输入姿势影响。以上指标均不能给出结合自由能、酯化初速、热稳定性或合成/水解比。

## 结构复核

工作台按设计序列的 SHA-256 查找本机已完成的 AF2 任务。选中相符的 rank_001 模型后可一键复核，记录预测编号、pTM、平均 pLDDT 和可用的三联体局部 pLDDT。也可上传候选对应的模型或实验结构，或使用命令行。将下例的批次、候选标识、坐标及链替换为实际值。原子坐标需与本条候选序列精确一致，否则程序给出原因。设计是亲本的少量替换，因此采用成对比对，保留高一致性结果。

~~~sh
python -m bur01_structmine.cli design-review --run results/designs/bur01_example --candidate AD_ID --structure candidate.pdb --chain A --origin predicted
~~~

坐标来源可标为 experimental、predicted、template_model 或 unspecified。模板建模保留亲本骨架时，极小的主链 RMSD 反映建模方式。对折叠或局部构象的独立判断需要另行获得的结构，对辛酸乙酯合成及相应水解的判断需要同条件酶学测定。
