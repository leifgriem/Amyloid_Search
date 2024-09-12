import os
import shutil
import argparse
from tqdm import tqdm

# Set up argument parsing
parser = argparse.ArgumentParser(description='Copy rank_001 PDB files from one directory to another.')
parser.add_argument('--alphafold_raw_output', type=str, required=True, help='Directory containing raw AlphaFold output.')
parser.add_argument('--output_dir', type=str, required=True, help='Directory to copy the rank_001 PDB files to.')
args = parser.parse_args()

# Access the arguments
alphafold_raw_output = args.alphafold_raw_output
output_dir = args.output_dir

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Get a list of all files in the specified directory
files_to_copy = [filename for filename in os.listdir(alphafold_raw_output) if 'rank_001' in filename and filename.endswith('.pdb')]

# Iterate over all files and copy them, using tqdm for progress tracking
for filename in tqdm(files_to_copy, desc="Copying rank_001 PDB files"):
    # Construct full file paths
    src_file = os.path.join(alphafold_raw_output, filename)
    dst_file = os.path.join(output_dir, filename)
    
    # Copy the file
    shutil.copy(src_file, dst_file)


