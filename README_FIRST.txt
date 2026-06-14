Bur01-StructMine V1 PyCharm Project

Please open the whole folder "Bur01_StructMine_PyCharm" in PyCharm.
All file names are English and all text files are UTF-8 encoded.

Fast test without Foldseek:
1. Open this folder in PyCharm.
2. Create Python interpreter: Python 3.10 or 3.11.
3. In PyCharm terminal, run:
   pip install -r requirements.txt
   python run_demo_no_foldseek.py
4. If results/demo_query_sequence.fasta is generated, the project is installed correctly.

Full run with Foldseek:
1. Use Linux or Windows WSL2.
2. Install Foldseek:
   conda install -c conda-forge -c bioconda foldseek -y
3. Download database:
   bash scripts/download_pdb100_database.sh
4. Run:
   python run_pipeline.py --config config.yaml

Demo with the included example Foldseek TSV:
python run_pipeline.py --config config.yaml --skip-foldseek --input-tsv data/example_foldseek_result.tsv
