import pandas as pd
import numpy as np
import os
import sys

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("Running convert_bim_idto_chr_pos.py. Usage: python script.py bfile_prefix out")
        print("Example: python3 convert_bim_idto_chr_pos.py ...")
        sys.exit(1)
    bim_file = sys.argv[1]
    input_dir = sys.argv[2]
    
    file_name = f"{bim_file}.bim"
    file_path = os.path.join(input_dir, file_name)
    if not os.path.exists(file_path):
        print("bim file path not found, snp id convertation halted")

    
    # Modifying snp IDs in bim file (rsIDs into chr:pos format)
    bim_df = pd.read_csv(file_path, sep="\t", header=None, names=["CHR", "RSID", "CM", "POS", "A1", "A2"])
    bim_df["RSID"] = bim_df["CHR"].astype(str) + ":" + bim_df["POS"].astype(str)
    modified_bim_file = f"{bim_file.split('.')[0]}_modified.bim"
    bim_df.to_csv(os.path.join(input_dir, f"{file_name}_modified.bim"), sep="\t", header=False, index=False)
    print(f"Modified BIM file saved to f'{file_name}_modified.bim'")
    
