# 0.5.0.dev0 本地验证 / Local Validation

日期 / Date：2026-09-21。该记录只覆盖本轮本地改造，不代表完整 NR、真实酶学性能或正式发布验收。

## 自动测试 / Automated Tests

- 96 项通过 / 96 passed，15.77 s；记录：[pytest.xml](validation/v0.5/pytest.xml)。所有测试使用独立临时目录，无正式 results 污染。
- Python 3.11.9；`pip check`：No broken requirements found；源代码、网页和测试编译通过。
- 最终文件检查 / Final file check：119个必需文件全部存在；已安装包元数据与项目配置均为0.5.0.dev0。
- 覆盖评分公式、Ser中心模体、缺失侧链、历史字段兼容、独立重新评分/重新注释、导入错误、重复测量、条件/单位分组、冲突候选、真实R²、相关性限制、Top-k、模型/链/插入码和只读演示。
- 原始回归中4项Windows批处理调用在含空格路径下失败；修复cmd参数引号后恢复。未删除或放宽这些测试。
- 本地可编辑包版本已更新为0.5.0.dev0。初次非隔离构建缺少bdist_wheel，改用隔离构建后成功；临时构建及缓存位于D盘。

## 真实 Foldseek / Real Foldseek

命令 / Command：

```powershell
python -m bur01_structmine.cli benchmark --engine foldseek --no-prepare --output-root results/validation_v05/benchmark
```

运行编号：`20260921_020257_509329_Bur01_benchmark_foldseek`。
Foldseek：`7f4b94644d6b210af4530e04a9fdc9e4bb64ea66`，通过Ubuntu-Foldseek实际运行，未使用Python替代正式检索。

| 检查 / Check | 结果 / Result |
|---|---|
| 8HFW自比对 qTM / tTM | 1.000 / 1.000 |
| 自比对RMSD | 0.000 Å |
| 所有阳性高于两个阴性 | 通过；原始alntmscore最低阳性0.7842，最高阴性0.2236 |
| P00520默认结构过滤 | 拒绝 / rejected |
| 原始近似alntmscore | 自比对为1.004；保留原值，不误当成归一化qTM/tTM |

详细清单与原始输出保留在对应运行目录。公开演示索引：[demo_manifest.json](../data/demo/demo_manifest.json)。当前源码重算演示评分与已保存值一致，最大分差0。

## Swiss-Prot 重新评分 / Rescoring

- 原运行 / Source：`20260711_AF2_8HFW_SwissProt_final`，12条候选。
- 新记录 / New record：`results/validation_v05/rescores/20260921_020553_201328_rescore/`。
- 原表SHA-256：`4b3da4b1b25825da6fee5e9c7305d8bdbdaeedcb54b88938935c8501b2a03179`。
- 11条因旧版非定位模体加分取消而降低5分，1条不变。未增加网络下载或推测缺失坐标；原预测运行未记录催化参考映射的候选明确保留此缺项。
- 原表及清单未改写，新清单记录源运行、输入哈希、版本和分差。本轮没有重跑完整Swiss-Prot搜索，重新评分使用已有真实搜索记录。

## 实验统计 / Experimental Statistics

模拟示例与真实数据分目录保存。5个样品，每个2次重复：均值10/20/30/40/50，样本标准差sqrt(2)。模拟预测12/18/33/38/49。

MAE=2，RMSE=sqrt(4.4)，R²=0.978；WASE排序Spearman=0.9，成功阈值≥30、k=3时WASE命中率2/3。另以y=[1,2,3]、p=[11,12,13]验证R²=−149而Pearson=1，排除了相关系数平方替代R²的错误。

命令行真实读写结果：`results/validation_v05/experiment_cli/`。这些数值全部为模拟演示，不是酶学验证结果。

## 网页与停止检查 / UI and Shutdown

- 通过 / PASS：1440×900桌面和390×844手机视口，中文与英文各检查5个演示步骤，共20个状态；均无页面横向溢出，键盘焦点落在可操作控件。原始记录：[ui_report.json](validation/v0.5/ui_report.json)。
- 三维画布像素非空，演示画布最低像素方差2176.6894（检测阈值3）；鼠标拖动前后的截图哈希发生变化，确认结构可以旋转。标签调整为选择Ser/Asp/His后显示单个位点，避免多个标签重叠。
- 工作台完成桌面/手机×中文/英文共4组检查：实际上传模拟CSV并导出Excel；用空格键切换Top-k；散点和排序图像素非空。上传公开1XKL结构，选择已有STH配体，并将A:GLU3标记为“关注”；四组JSON导出均包含1个位点标记及坐标存在状态。
- 断网展示采用新浏览器上下文阻断所有非本机地址请求，记录到的外部请求为0；后端隔离测试同时禁止演示调用网络和计算子进程。这是受控网络隔离测试，未声称物理断开网线。三维脚本、结构和结果均从本地读取。
- 工作台检查使用独立固定基准记录和测试模式诊断，不代表新增GPU预测或完整NR验收。实验示例始终标为模拟数据。
- 测试结束后已停止演示和工作台服务，终止本轮启动的Ubuntu-Foldseek；再次检查未发现项目Python服务、devtunnel或运行中的WSL。未增加自启动、定时任务或公网通道。演示按需通过`start_demo.bat`启动，通过`stop_demo.bat`停止。

English summary: all 20 demo states and four bilingual desktop/mobile workbench combinations passed. Checks covered keyboard interaction, CSV import/Excel export, site-mark JSON export, nonblank charts and 3D canvases, and structure rotation. Non-local requests were blocked in fresh browser contexts with zero external requests observed; this was controlled network isolation, not a physical disconnection. Test services and the WSL instance used for the benchmark were stopped. No project Python service, devtunnel process or running WSL distribution remained at the final check. These checks do not constitute GPU-prediction, experimental-function or full-NR acceptance.

截图示例 / Example captures：[桌面演示](validation/v0.5/demo_desktop_zh_3.png)、[手机演示](validation/v0.5/demo_mobile_zh_3.png)、[中文实验评估](validation/v0.5/workbench_desktop_zh_experiments.png)、[英文手机位点查看](validation/v0.5/workbench_mobile_en_sites.png)。

## 保留与未完成项 / Preserved and Pending

- 保留已有结果、数据库、Obsidian库、交付包和进入本轮前的未提交改动（含benchmark下载清单与patent目录）。
- 完整NR服务器导入验收仍未完成，状态未因界面改造变更。
- 未新增序列生成、酶性能预测模型或真实实验数据；未重新进行GPU AlphaFold预测。
- 本轮不推送GitHub，不合并PR，不创建Release，也不部署公网服务。
