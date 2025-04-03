import pandas as pd
import numpy as np
import sys
import os

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("Usage: python script.py <input_dir> <file_prefix>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    file_prefix = sys.argv[2]
    
    file_path = os.path.join(input_dir, f"{file_prefix}.0.3.profile")
    if not os.path.exists(file_path):
        print("profile file path not found, z-normalization halted")
    
    profile = pd.read_csv(file_path, delim_whitespace=True)
    profile['SCORE'] = (profile['SCORE'] - profile['SCORE'].mean()) / profile['SCORE'].std()
    profile.to_csv(f"{file_prefix}_0.3.znorm.profile", sep="\t", index=False)
