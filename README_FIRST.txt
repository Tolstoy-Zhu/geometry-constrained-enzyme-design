AquaEsterNet-StructMine v0.3 快速入口
=====================================

1. 打开网页界面
   双击 start_ui.bat
   固定地址：http://127.0.0.1:8501

   首次启用“登录后自动启动”：双击 install_ui_autostart.bat
   查看状态：双击 ui_status.bat
   停止网页：双击 stop_ui.bat

2. 在线展示给其他电脑
   首次双击 install_online_ui.bat
   分享文件：.cache\remote_ui\SHARE_WITH_EXPERTS.txt
   停止公网访问：双击 stop_online_ui.bat

3. 环境诊断
   python start_here.py --action doctor

4. 固定基准的正式 Foldseek 验收
   python start_here.py --action benchmark

5. 完整 AlphaFold/Swiss-Prot 搜索
   python start_here.py --action production-run

6. 单体 AlphaFold2 预测
   python -m bur01_structmine.cli predict --fasta query.fasta --preset standard --confirm-public-msa

7. 预测后直接结构挖掘
   python -m bur01_structmine.cli predict-and-mine --fasta query.fasta --preset standard --confirm-public-msa --mining-config configs/alphafold_swissprot.yaml

8. 运行测试
   python start_here.py --action validate

9. 检查项目文件
   python start_here.py --action check-files

查询结构：PDB 8HFW chain A，残基 12-294，催化残基 Ser111/Asp241/His274。

重要说明：
- target-run 现在等同于真实 Foldseek benchmark。
- P00520 是远缘负对照，不是阳性候选。
- TM-align 仅用于小集合预览，不能替代 Foldseek 最终验收。
- WASE-screen 是透明预筛分，不预测 S/H ratio，也不等同于正式 WASE。
- 每次运行写入 results/runs/<run_id>/，失败和取消均有清单记录。

本机最终验收：43 项测试通过，31 项诊断通过且无警告；本地网页与在线 HTTPS 通道的崩溃恢复通过；8HFW 和 Q9FVW3 两套 5 模型 AF2 标准预测通过；预测 8HFW 的正式 Foldseek 基准通过；完整 Swiss-Prot v6 2025_03 搜索得到 999 个原始命中、12 个过滤候选和 12 个有效注释；桌面与手机网页验收通过。

详细中文说明：README_CN.md
模型原理：MODEL_GUIDE_CN.md
验证证据：VALIDATION_REPORT.md
