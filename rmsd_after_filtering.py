# rmsd_after_filtering.py
# Run an additional RMSD calculation for refined list of Alphafold structures generated by Alignment_MDAnalysis.py

import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Run RMSD calculations for refined list of AlphaFold structures.")

# Add arguments
parser.add_argument('--alphafold_dir', type=str, required=True, help='Directory containing AlphaFold PDB files.')
parser.add_argument('--amyloid_dir', type=str, required=True, help='Directory containing amyloid PDB files.')
parser.add_argument('--csv_dir', type=str, required=True, help='Directory to save CSV files.')
parser.add_argument('--img_dir', type=str, required=True, help='Directory to save images.')
parser.add_argument('--protein', type=str, required=True, help='Name of the protein being analyzed.')
parser.add_argument('--segment_length', type=int, required=True, help='Segment length.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
alphafold_dir = args.alphafold_dir
amyloid_dir = args.amyloid_dir
csv_dir = args.csv_dir
img_dir = args.img_dir
protein = args.protein
segment_length = args.segment_length

# Ensure output directories exist
os.makedirs(csv_dir, exist_ok=True)
os.makedirs(img_dir, exist_ok=True)

# Dynamically find n_range by parsing through files in alphafold_dir
n_range = 0
# Use glob to get a list of all PDB files in the AlphaFold directory
alphafold_files = glob.glob(os.path.join(alphafold_dir, "*.pdb"))
# Find largest resid value to use as n_range
for alphafold_file in alphafold_files:
    # Get the base name of the file
    base_name = os.path.basename(alphafold_file)
    
    # Split the base name by '_'
    parts = base_name.split('_')
    
    # Search for 'digits-digits' pattern using regex
    for part in parts:
        match = re.match(r"(\d+)-(\d+)", part)
        if match:
            # Extract the first digit
            first_digit = int(match.group(1))
            # Update n_range if the first digit is greater than current n_range
            if first_digit > n_range:
                n_range = first_digit

# Use glob to get a list of all PDB files in the amyloid directory
amyloid_files = glob.glob(os.path.join(amyloid_dir, "*.pdb"))

# Initialize a dataframe to store RMSD values with all values set to np.inf
df = pd.DataFrame(np.inf, index=[f"{n}-{n + segment_length - 1}" for n in range(1, n_range + 1)],
                  columns=[os.path.splitext(os.path.basename(amyloid_file))[0] for amyloid_file in amyloid_files])

# Loop through each AlphaFold PDB file with tqdm progress bar
for alphafold_file in tqdm(alphafold_files, desc="Processing AlphaFold PDB files"):

    # Extract just the base name of the file
    base_name = os.path.basename(alphafold_file)
    
    # Split the base name to extract the amyloid PDB name and residue range
    file_parts = base_name.split('_')
    amyloid_name = file_parts[0]  # The amyloid PDB name
    residue_range = file_parts[1]
    first_resid, last_resid = map(int, residue_range.split('-'))

    # Construct the path for the corresponding amyloid PDB file
    amyloid_file = os.path.join(amyloid_dir, f"{amyloid_name}.pdb")

    # Ensure the corresponding amyloid PDB file exists
    if not os.path.exists(amyloid_file):
        print(f"Corresponding amyloid PDB file {amyloid_file} not found. Skipping...")
        continue

    # Load universes
    alphafold_universe = mda.Universe(alphafold_file)
    amyloid_universe = mda.Universe(amyloid_file)

    # Define selection string for the AlphaFold PDB
    alphafold_selection = f"chainID A and resid {first_resid}:{last_resid} and (name N CA C O CB)"
    alphafold_selection_without_cb = f"chainID A and resid {first_resid}:{last_resid} and (name N CA C O)"
    alphafold_atoms = alphafold_universe.select_atoms(alphafold_selection)

    # Define selection string for the amyloid PDB (backbone atoms + CB)
    amyloid_selection = f"chainID A and resid {first_resid}:{last_resid} and (name N CA C O CB)"
    amyloid_selection_without_cb = f"chainID A and resid {first_resid}:{last_resid} and (name N CA C O)"
    amyloid_atoms = amyloid_universe.select_atoms(amyloid_selection)

    # Check if the selections have different numbers of atoms
    if len(amyloid_atoms) != len(alphafold_atoms):
        # Re-select without CB atoms
        alphafold_atoms = alphafold_universe.select_atoms(alphafold_selection_without_cb)
        amyloid_atoms = amyloid_universe.select_atoms(amyloid_selection_without_cb)

        # Skip if the new selections still don't match
        if len(amyloid_atoms) != len(alphafold_atoms):
            print(f"Selections do not match for {base_name}. Skipping RMSD calculation.")
            continue

    # Calculate RMSD
    rmsd_value = rms.rmsd(alphafold_atoms.positions, amyloid_atoms.positions, center=True, superposition=True)
    # Convert from Angstroms to Nanometers
    rmsd_value = rmsd_value / 10

    # Define row and column names explicitly
    row_name = f"{first_resid}-{last_resid}"
    pdb_name = amyloid_name  # Use the extracted amyloid PDB name

    # Store the RMSD value in the dataframe
    df.loc[row_name, pdb_name] = rmsd_value

# Replace NaN values with np.inf in the dataframe
df = df.fillna(np.inf)
#print(df)
#%%
# Save dataframe as csv in the specified directory
csv_output_path = os.path.join(csv_dir, f'RMSD_values_final_{protein}_{segment_length}aa.csv')
df.to_csv(csv_output_path)

# Generate plots to analyze data
#%%
# Generate a heatmap of the RMSD values
plt.figure(figsize=(12, 8))
sns.heatmap(df, annot=False, cmap='viridis', cbar=True, vmin=0.2, vmax=0.4)
plt.title(f'Filtered RMSD Heatmap {protein} {segment_length}aa')
plt.xlabel('Amyloid PDB')
plt.ylabel(f'Residue Range (n:n+{segment_length-1})')

# Save heatmap as an image in specified directory
heatmap_output_path = os.path.join(img_dir, f"{protein}_{segment_length}aa_Filtered_RMSD_Heatmap.png")
plt.savefig(heatmap_output_path)
#plt.show()

#%%
# Generate a histogram of RMSD values (ignoring np.inf)
plt.figure(figsize=(8, 6))
rmsd_values = df.replace(np.inf, np.nan).values.flatten()
rmsd_values = rmsd_values[~np.isnan(rmsd_values)]  # Remove NaNs
plt.hist(rmsd_values, bins=30, color='blue', alpha=0.7)
plt.title(f'Filtered RMSD Histogram {protein} {segment_length}aa')
plt.xlabel('RMSD (nm)')
plt.ylabel('Frequency')

# Save histogram as an image in specified directory
histogram_output_path = os.path.join(img_dir, f"{protein}_{segment_length}aa_Filtered_RMSD_Histogram.png")
plt.savefig(histogram_output_path)
#plt.show()
