AquaEsterNet-StructMine v0.4 快速入口
=====================================

1. 临时打开本地网页
   双击 start_ui.bat
   地址：http://127.0.0.1:8501
   使用完毕后双击 stop_ui.bat。

   本轮 v0.4 不安装、不恢复后台自启动，也不自动打开公网隧道。

2. 环境诊断
   python start_here.py --action doctor

3. 结构挖掘
   python start_here.py --action benchmark
   python start_here.py --action production-run

4. AlphaFold2 单体预测
   python -m bur01_structmine.cli predict --fasta query.fasta --preset standard --confirm-public-msa

5. 构建真实 NCBI NR 本机命中子集
   需要联系邮箱，并逐次确认 Bur01 序列会发送至 NCBI：
   python -m bur01_structmine.cli nr-demo-build --email name@example.org --confirm-ncbi

6. 从结构运行扩增 NR
   python -m bur01_structmine.cli nr-expand --source-run <run_id> --config configs/nr_local_demo.yaml

7. 完整 NR 服务器/HPC
   python -m bur01_structmine.cli nr-export-job --source-run <run_id> --config configs/nr_full_server.yaml
   python -m bur01_structmine.cli nr-import-job --bundle <nr_run_id>_server_result.zip

8. NR 新颖性审计
   python -m bur01_structmine.cli nr-novelty --fasta designed.fasta --backend ncbi_remote --email name@example.org --confirm-ncbi

9. 测试与文件完整性
   python start_here.py --action validate
   python start_here.py --action check-files

重要边界：
- P00520 是远缘负对照，不是阳性候选。
- NR 序列命中不是结构相似性证据。
- 未完成 AF2 和 Foldseek/TM-align 复核的 NR 候选保持 WASE not evaluated。
- WASE-screen 不预测 S/H ratio，也不等同于正式 WASE。
- 完整 NR 结果导回并验收前，不合并 v0.4 分支、不打标签、不创建 Release。

详细中文说明：README_CN.md
NR 操作说明：docs/NR_WORKFLOW_CN.md
模型原理：MODEL_GUIDE_CN.md
验证状态：VALIDATION_REPORT.md
