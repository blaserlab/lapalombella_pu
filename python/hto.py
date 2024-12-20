import scanpy as sc
import pandas as pd
from pathlib import Path

# input_folder = Path("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/fastq_pipestance/")
input_folder = Path("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/fastq_pipestance_batch2/")
output_folder = Path("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/staff/single_cell/lapalombella_citeseq_dec2024/hto_lookup/")

files = [
    # (input_folder / "PZ-3606_WT-pool1_CITE/Hashtag_results/adataFinal_PZ-3606_WT-pool1.h5ad", output_folder / "WT_pool1.csv"),
    # (input_folder / "PZ-3625_WT-pool2_CITE/Hashtag_results/adataFinal_PZ-3625_WT-pool2.h5ad", output_folder / "WT_pool2.csv"),
    # (input_folder / "PZ-3630_COMUTANT-pool2_CITE/Hashtag_results/adataFinal_PZ-3630_COMUTANT-pool2.h5ad", output_folder / "comutant_pool2.csv"),
    # (input_folder / "PZ-3630_TET2-pool2_CITE/Hashtag_results/adataFinal_PZ-3630_TET2-pool2.h5ad", output_folder / "tet2_pool2.csv"),
    # (input_folder / "PZ-3630_TP53-pool2_CITE/Hashtag_results/adataFinal_PZ-3630_TP53-pool2.h5ad", output_folder / "tp53_pool2.csv")
    (input_folder / "PZ-3616_COMUTANT-pool1_CITE/Hashtag_0.1.3_results/outputs/PZ-3616_COMUTANT-pool1.h5ad", output_folder / "comutant_pool1.csv"),
    (input_folder / "PZ-3616_TET2-pool1_CITE/Hashtag_0.1.3_results/outputs/PZ-3616_TET2-pool1.h5ad", output_folder / "tet2_pool1.csv"),
    (input_folder / "PZ-3616_TP53-pool1_CITE/Hashtag_0.1.3_results/adataFinal_PZ-3616_TP53-pool1.h5ad", output_folder / "tp53_pool1.csv"),
    (input_folder / "PZ-3641_COMUTANT-pool3_CITE/Hashtag_0.1.3_results/outputs/PZ-3641_COMUTANT-pool3.h5ad", output_folder / "comutant_pool3.csv"),
    (input_folder / "PZ-3641_COMUTANT-pool4_CITE/Hashtag_0.1.3_results/outputs/PZ-3641_COMUTANT-pool4.h5ad", output_folder / "comutant_pool4.csv")
    ]

for input_file, output_file in files:
    hto = sc.read_h5ad(input_file)
    hto_obs = pd.DataFrame(hto.obs)
    hto_obs.to_csv(output_file, index = False)
