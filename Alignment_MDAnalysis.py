import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align
import os
import warnings
import pandas as pd
import numpy as np
import re
import glob
from tqdm import tqdm

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Align structures and calculate interchain contacts and RMSDs.")

# Add arguments for input file paths, directories, and other parameters
parser.add_argument('--interchain_contacts', type=str, required=True, help='Path to the interchain contacts CSV file.')
parser.add_argument('--RMSDs', type=str, required=True, help='Path to the RMSDs CSV file.')
parser.add_argument('--alphafold_dir', type=str, required=True, help='Directory containing AlphaFold PDB files.')
parser.add_argument('--amyloid_dir', type=str, required=True, help='Directory containing amyloid PDB files.')
parser.add_argument('--output_dir', type=str, required=True, help='Directory to save aligned structures.')
parser.add_argument('--protein', type=str, required=True, help='Name of the protein being analyzed.')
parser.add_argument('--segment_length', type=int, required=True, help='Segment length.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
interchain_contacts = args.interchain_contacts
RMSDs = args.RMSDs
alphafold_dir = args.alphafold_dir
amyloid_dir = args.amyloid_dir
output_dir = args.output_dir
protein = args.protein
segment_length = args.segment_length

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Interchain contacts threshold
interchain_threshold = 0.7
# 175 = 
RMSD_threshold = 0.35

# Suppress specific MDAnalysis warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Unit cell dimensions not found.")
warnings.filterwarnings("ignore", category=UserWarning, message="Found no information for attr: 'formalcharges'")

def alignment(alphafold_pdb, amyloid_pdb, output_dir):
    # Select MDanalyis universes for ref and mobile
    amyloid_universe = mda.Universe(amyloid_pdb)
    alphafold_universe = mda.Universe(alphafold_pdb)

    # Extract the Amyloid name without .pdb attached
    amyloid_name = os.path.splitext(os.path.basename(amyloid_pdb))[0]

    # Extract the numbers in the alphafold_pdb filename
    alphafold_name = os.path.basename(alphafold_pdb)
    parts = alphafold_name.split('_')

    # Initialize firstnum and lastnum
    firstnum, lastnum = None, None

    # Apply regex to each part of the split filename
    for part in parts:
        match = re.match(r'(\d+)-(\d+)', part)
        if match:
            firstnum, lastnum = map(int, match.groups())
            break

    # Select atoms firstnum to lastnum for ref universe
    amyloid_selection = f"chainid A and resid {firstnum}:{lastnum} and name N CA C"
    # Generalize to {2*segment_length+21}:{3*segment_length+20}
    alphafold_selection = f"chainid A and resid {2*segment_length+21}:{3*segment_length+20} and name N CA C"
    # Assign selection to mobile_atoms and ref_atoms
    amyloid_atoms = amyloid_universe.select_atoms(f"{amyloid_selection}")
    alphafold_atoms = alphafold_universe.select_atoms(f"{alphafold_selection}")

    # If len(amyloid_atoms) == len(alphafold_atoms), then align the structures
    if len(amyloid_atoms) == len(alphafold_atoms):

        # Perform MDanalysis align.alignto command
        align.alignto(alphafold_universe, amyloid_universe, select=(alphafold_selection, amyloid_selection), weights="mass")

        # Select atoms to keep in output_universe
        first_resid, last_resid = segment_length+11, 4*segment_length+30
        atoms_to_keep = f"chainid A and resid {first_resid}:{last_resid}"
        selected_atoms = alphafold_universe.select_atoms(atoms_to_keep)

        # Name the aligned structure "{amyloid_name}_{firstnum}-{lastnum}_aligned.pdb"
        aligned_filename = f"{amyloid_name}_{firstnum}-{lastnum}_aligned.pdb"
        # Save the aligned structure to output directory
        aligned_path = os.path.join(output_dir, aligned_filename)

        # Write aligned structure with selected residues to output directory
        with mda.Writer(aligned_path, selected_atoms.n_atoms) as W:
            W.write(selected_atoms)

# Read interchain_connects and RMSDs as dataframes
interchain_df = pd.read_csv(interchain_contacts, index_col=0)
RMSD_df = pd.read_csv(RMSDs, index_col=0)

# Step 1: Ensure indices are consistent between both DataFrames
RMSD_df.index = RMSD_df.index.astype(str)
interchain_df.index = interchain_df.index.astype(str)

# Trim the larger DataFrame to match the smaller one by row count
if len(RMSD_df) > len(interchain_df):
    RMSD_df = RMSD_df.iloc[:len(interchain_df)]
elif len(interchain_df) > len(RMSD_df):
    interchain_df = interchain_df.iloc[:len(RMSD_df)]

# Step 3: Apply the interchain mask to filter the RMSD DataFrame
interchain_mask = interchain_df['Interchain Contacts'] > interchain_threshold

# Convert DataFrames to NumPy arrays
RMSD_array = RMSD_df.to_numpy()
interchain_mask_array = interchain_mask.to_numpy()

# Apply the mask to filter rows in RMSD_df
RMSD_filtered_array = RMSD_array[interchain_mask_array]

# Convert the filtered array back to a DataFrame, keeping the original column names
RMSD_df_filtered = pd.DataFrame(RMSD_filtered_array, columns=RMSD_df.columns, index=RMSD_df.index[interchain_mask_array])

# Save results as a list of tuples
alphafold_and_amyloid = [(row, col) for (row, col), value in RMSD_df_filtered.stack().items() if value < RMSD_threshold]

# Sort alphafold_and_amyloid alphabetically by the second tuple value for every entry
alphafold_and_amyloid.sort(key=lambda x: x[1])

# Iterate through each tuple (resid, pdb) in resid_and_pdb
for alphafold, amyloid in tqdm(alphafold_and_amyloid, desc="Processing pairs"):
    # Replace underscores with hyphens for the alphafold filenames
    alphafold = alphafold.replace('_', '-')

    # Iterate through amyloid_dir that contains {amyloid}
    found_any_files = False
    for amyloid_file in glob.glob(os.path.join(amyloid_dir, f'*{amyloid}*')):
        # Iterate through alphafold_dir that contains {alphafold}
        for alphafold_file in glob.glob(os.path.join(alphafold_dir, f'*{alphafold}*')):
            found_any_files = True
            # Align the structures alphafold_file, amyloid_file, output_dir
            alignment(alphafold_file, amyloid_file, output_dir)
    if not found_any_files:
        print(f"No matching files found for pair: {alphafold}, {amyloid}")

