# Initial_RMSD.py
# Calculate RMSD between AlphaFold predictions and reference PDB structures using MDAnalysis. Store results as a dataframe and optionally generate a heatmap.
#%% 

import argparse
import os
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import pandas as pd
import numpy as np
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Calculate RMSD between AlphaFold predictions and reference PDB structures using MDAnalysis.")

# Add arguments
parser.add_argument('--alphafold_dir', type=str, required=True, help='Directory containing AlphaFold simulations.')
parser.add_argument('--ref_dir', type=str, required=True, help='Directory containing reference PDB files.')
parser.add_argument('--csv_dir', type=str, required=True, help='Directory to save the RMSD dataframe as .csv.')
parser.add_argument('--img_dir', type=str, required=True, help='Directory to save heatmap of results.')
parser.add_argument('--protein', type=str, required=True, help='Name of the protein being analyzed.')
parser.add_argument('--segment_length', type=int, required=True, help='Segment length.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
alphafold_dir = args.alphafold_dir
ref_dir = args.ref_dir
csv_dir = args.csv_dir
img_dir = args.img_dir
protein = args.protein
segment_length = int(args.segment_length)

# Create the output directories if they don't exist
os.makedirs(csv_dir, exist_ok=True)
os.makedirs(img_dir, exist_ok=True)

def find_chain_residues(pdb_path):
    # Load the reference structure
    ref = mda.Universe(pdb_path)
    # Select residues from chain A using ChainID
    chain_a_residues = ref.select_atoms("chainid A").residues
    resid_ids = chain_a_residues.resids
    
    if len(resid_ids) == 0:
        return None, None  # No residues found in chain A

    # Find the first and last residue IDs
    first_resid = resid_ids[0]
    last_resid = resid_ids[-1]
    
    return first_resid, last_resid

def extract_range(filename):
    '''
    filename: str, the name of the pdb file
    '''
    # Extract the ##-## part from the filename and return the first number
    # Regular expression to find patterns like "digits-digits" (e.g., "32-46")
    match = re.search(r'(\d+)-(\d+)', filename)
    if match:
        start_num = int(match.group(1))  # Extract the first number
        return start_num
    else:
        return -1  # Return a default value if no match is found

# Create alphafold_dict and pdb_dict
alphafold_files = [f for f in os.listdir(alphafold_dir) if "rank_001" in f and f.endswith(".pdb")]
sorted_alphafold_files = sorted(alphafold_files, key=extract_range)
alphafold_dict = {extract_range(f): os.path.join(alphafold_dir, f) for f in sorted_alphafold_files}

ref_files = [f for f in os.listdir(ref_dir) if f.endswith(".pdb")]
pdb_dict = {i: os.path.join(ref_dir, f) for i, f in enumerate(ref_files, start=1)}

# Initialize a dataframe to store RMSD values
height = max(alphafold_dict.keys()) + 1
length = len(pdb_dict)+1
df = pd.DataFrame(np.full((height, length), np.inf))

column_names = ["PDB ID"] + [os.path.splitext(os.path.basename(path))[0] for path in pdb_dict.values()]
df.columns = column_names

df["PDB ID"] = [f"{i}_{i+segment_length-1}" for i in range(1, height + 1)]
pdb_name_dict = {key: os.path.splitext(os.path.basename(path))[0] for key, path in pdb_dict.items()}

# Start RMSD calculations
for key, ref_path in tqdm(pdb_dict.items(), desc="Processing reference files"):
    ref_universe = mda.Universe(ref_path)
    first_resid, last_resid = find_chain_residues(ref_path)

    if first_resid is None or last_resid is None:
        continue  # Skip if no residues found

    for i in range(first_resid, last_resid - segment_length + 1):
        if i not in alphafold_dict:
            continue

        mobile_path = alphafold_dict[i]
        mobile_universe = mda.Universe(mobile_path)

        mobile_selection = f"chainID A and resid {2*segment_length+21}:{3*segment_length+20} and name CA"
        ref_selection = f"chainID A and resid {i}:{i + segment_length - 1} and name CA"

        ref_atoms = ref_universe.select_atoms(f"{ref_selection}")
        mobile_atoms = mobile_universe.select_atoms(f"{mobile_selection}")

        # Skip RMSD calc if lengths don't match
        if len(ref_atoms) != len(mobile_atoms):
            ref_atoms = ref_universe.select_atoms(f"chainID A and resid {i}:{i + segment_length - 1} and name CA")
            mobile_atoms = mobile_universe.select_atoms(f"chainID A and resid {2*segment_length+21}:{3*segment_length+20} and name CA")
            continue

        rmsd_value = rms.rmsd(mobile_atoms.positions, ref_atoms.positions, center=True, superposition=True)
        rmsd_value = rmsd_value / 10
        df.loc[i, pdb_name_dict[key]] = rmsd_value

# Save the dataframe to a csv file in the user specified directory
df.to_csv(os.path.join(csv_dir, f"initial_RMSD_{protein}_{segment_length}aa.csv"), index=False)

# Generate heatmap
plt.figure(figsize=(10, 8))
mask = df.iloc[:, 1:] <= 0
sns.heatmap(df.iloc[:, 1:], mask=mask, cmap="viridis", cbar_kws={"label": "RMSD (nm)"}, annot=False, vmin=0.1, vmax=0.3)
plt.xlabel('Reference PDBs')
plt.ylabel(f'Alphafold Amyloid Predictions (n:n+{segment_length-1})')
plt.savefig(os.path.join(img_dir, f"Initial_RMSD_Heatmap_{segment_length}aa.png"))
#plt.show()

# Generate histogram
rmsd_values = df.iloc[:, 1:].to_numpy().flatten()
rmsd_values = rmsd_values[np.isfinite(rmsd_values)]
plt.figure(figsize=(8, 6))
plt.hist(rmsd_values, bins=30, color='blue', alpha=0.7, edgecolor='black')
plt.title(f'Histogram of RMSD Values {segment_length}aa')
plt.xlabel('RMSD (nm)')
plt.ylabel('Frequency')
plt.savefig(os.path.join(img_dir, f"Initial_RMSD_Histogram_{segment_length}aa.png"))
#plt.show()

# #%%
# import os
# import MDAnalysis as mda
# from MDAnalysis.analysis import rms
# import pandas as pd
# import numpy as np
# from tqdm import tqdm
# import seaborn as sns
# import matplotlib.pyplot as plt
# import re

# # Define variables directly in the script
# protein = "AB"  # Example value, change as necessary
# num = 13        # Example value, change as necessary
# segment_length = num

# # Directories using formatted strings with the variable `protein`
# alphafold_dir = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/{protein}/PDBs/{protein}_{num}aa_consecutive_filtered"
# ref_dir = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/{protein}/PDBs/{protein}_unique_amyloids"
# csv_dir = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/{protein}/CSVs"
# img_dir = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/{protein}/IMGs"

# # Create the output directories if they don't exist
# os.makedirs(csv_dir, exist_ok=True)
# os.makedirs(img_dir, exist_ok=True)

# def find_chain_residues(pdb_path):
#     ref = mda.Universe(pdb_path)
#     chain_a_residues = ref.select_atoms("chainid A").residues
#     resid_ids = chain_a_residues.resids
    
#     if len(resid_ids) == 0:
#         return None, None  # No residues found in chain A

#     first_resid = resid_ids[0]
#     last_resid = resid_ids[-1]
    
#     return first_resid, last_resid

# def extract_range(filename):
#     match = re.search(r'(\d+)-(\d+)', filename)
#     if match:
#         start_num = int(match.group(1))
#         return start_num
#     else:
#         return -1

# # Create alphafold_dict and pdb_dict
# alphafold_files = [f for f in os.listdir(alphafold_dir) if "rank_001" in f and f.endswith(".pdb")]
# sorted_alphafold_files = sorted(alphafold_files, key=extract_range)
# alphafold_dict = {extract_range(f): os.path.join(alphafold_dir, f) for f in sorted_alphafold_files}

# ref_files = [f for f in os.listdir(ref_dir) if f.endswith(".pdb")]
# pdb_dict = {i: os.path.join(ref_dir, f) for i, f in enumerate(ref_files, start=1)}

# # Initialize a dataframe to store RMSD values
# height = max(alphafold_dict.keys()) + 1
# length = len(pdb_dict) + 1
# df = pd.DataFrame(np.full((height, length), np.inf))

# column_names = ["PDB ID"] + [os.path.splitext(os.path.basename(path))[0] for path in pdb_dict.values()]
# df.columns = column_names
# df["PDB ID"] = [f"{i}_{i+segment_length-1}" for i in range(1, height + 1)]

# pdb_name_dict = {key: os.path.splitext(os.path.basename(path))[0] for key, path in pdb_dict.items()}

# # Start RMSD calculations
# for key, ref_path in tqdm(pdb_dict.items(), desc="Processing reference files"):
#     ref_universe = mda.Universe(ref_path)
#     first_resid, last_resid = find_chain_residues(ref_path)

#     if first_resid is None or last_resid is None:
#         continue  # Skip if no residues found

#     for i in range(first_resid, last_resid - segment_length + 1):
#         if i not in alphafold_dict:
#             continue

#         mobile_path = alphafold_dict[i]
#         mobile_universe = mda.Universe(mobile_path)

#         mobile_selection = f"chainID A and resid {2*segment_length+21}:{3*segment_length+20} and name CA"
#         ref_selection = f"chainID A and resid {i}:{i + segment_length - 1} and name CA"

#         ref_atoms = ref_universe.select_atoms(f"{ref_selection}")
#         mobile_atoms = mobile_universe.select_atoms(f"{mobile_selection}")

#         if len(ref_atoms) != len(mobile_atoms):
#             continue  # Skip if atom selections are not matching in length

#         rmsd_value = rms.rmsd(mobile_atoms.positions, ref_atoms.positions, center=True, superposition=True)
#         rmsd_value = rmsd_value / 10
#         df.loc[i, pdb_name_dict[key]] = rmsd_value

# # Save the dataframe to a csv file
# df.to_csv(os.path.join(csv_dir, f"initial_RMSD_{protein}_{segment_length}aa.csv"), index=False)

# # Generate heatmap
# plt.figure(figsize=(10, 8))
# mask = df.iloc[:, 1:] <= 0
# sns.heatmap(df.iloc[:, 1:], mask=mask, cmap="viridis", cbar_kws={"label": "RMSD (nm)"}, annot=False, vmin=0.1, vmax=0.3)
# plt.xlabel('Reference PDBs')
# plt.ylabel(f'Alphafold Amyloid Predictions (n:n+{segment_length-1})')
# plt.savefig(os.path.join(img_dir, f"Initial_RMSD_Heatmap_{segment_length}aa.png"))

# # Generate histogram
# rmsd_values = df.iloc[:, 1:].to_numpy().flatten()
# rmsd_values = rmsd_values[np.isfinite(rmsd_values)]
# plt.figure(figsize=(8, 6))
# plt.hist(rmsd_values, bins=30, color='blue', alpha=0.7, edgecolor='black')
# plt.title(f'Histogram of RMSD Values {segment_length}aa')
# plt.xlabel('RMSD (nm)')
# plt.ylabel('Frequency')
# plt.savefig(os.path.join(img_dir, f"Initial_RMSD_Histogram_{segment_length}aa.png"))



# # %%
