# Alignment_MDAnalysis.py
# Creates aligned PDB files with same shape as Alphafold predictions if RMSD < RMSD_threshold and Interchain Contacts > interchain_threshold
# %%
import MDAnalysis as mda
from MDAnalysis.analysis import align
import os
import warnings
import pandas as pd
import numpy as np

# Interchain contacts threshold
interchain_threshold = 175
RMSD_threshold = 0.35

# Suppress specific MDAnalysis warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Unit cell dimensions not found.")
warnings.filterwarnings("ignore", category=UserWarning, message="Found no information for attr: 'formalcharges'")

interchain_contacts = r"interchain_contacts_alpha_synuclein.csv"
RMSDs = r"RMSDs_Alpha_Synuclein.csv"

alphafold_dir = r"./PDBs/alpha_synuclein_15aa_consecutive_filtered"
amyloid_dir = r"./PDBs/alpha_synuclein_unique_amyloids"
output_dir = r"./PDBs/aligned_structures"

def alignment(alphafold_pdb, amyloid_pdb, output_dir):

    # Select MDanalyis universes for ref and mobile
    amyloid_universe = mda.Universe(amyloid_pdb)
    alphafold_universe = mda.Universe(alphafold_pdb)
    # Extract the Amyloid name without .pdb attached
    amyloid_name = os.path.splitext(os.path.basename(amyloid_pdb))[0]

    # Extract the numbers in the alphafold_pdb filename
    alphafold_name = os.path.basename(alphafold_pdb)
    alphafold_num = alphafold_name.split('_')[2]
    firstnum,lastnum = alphafold_num.split('-')
    firstnum,lastnum = int(firstnum),int(lastnum)

    # Select atoms firstnum to lastnum for ref universe
    amyloid_selection = f"chainid A and resid {firstnum}:{lastnum} and name N CA C"
    # Select atoms 51:65 for mobile universe
    alphafold_selection = "chainid A and resid 51:65 and name N CA C"
    # Assign selection to mobile_atoms and ref_atoms
    amyloid_atoms = amyloid_universe.select_atoms(f"{amyloid_selection}")
    alphafold_atoms = alphafold_universe.select_atoms(f"{alphafold_selection}")

    # If len(amyloid_atoms) == len(alphafold_atoms), then align the structures
    if len(amyloid_atoms) == len(alphafold_atoms):

        # Perform MDanalysis align.alignto command
        align.alignto(alphafold_universe, amyloid_universe, select=(alphafold_selection, amyloid_selection), weights = "mass")

        # Select atoms to keep in output_universe
        atoms_to_keep = "chainid A and resid 26:90"
        selected_atoms = alphafold_universe.select_atoms(atoms_to_keep)

        # Name the aligned structure "{amyloid_name}_{firstnum}-{lastnum}_aligned.pdb"
        aligned_filename = f"{amyloid_name}_{firstnum}-{lastnum}_aligned.pdb"
        # Save the aligned structure to output directory
        aligned_path = os.path.join(output_dir, aligned_filename)

        # Write aligned structure with selected residues to output directory
        with mda.Writer(aligned_path, selected_atoms.n_atoms) as W:
            W.write(selected_atoms)



#%%

# Read interchain_connects and RMSDs as dataframes
interchain_df = pd.read_csv(interchain_contacts,index_col=0)
RMSD_df = pd.read_csv(RMSDs,index_col=0)

#%%
# Step 1: Ensure indices are consistent between both DataFrames
RMSD_df.index = RMSD_df.index.astype(str)
interchain_df.index = interchain_df.index.astype(str)

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
# print(alphafold_and_amyloid)
#%%
import glob
# Iterate through each tuple (resid, pdb) in resid_and_pdb
for alphafold, amyloid in alphafold_and_amyloid:
    # Split resid into list of strings separated by _
    alphafold = alphafold.split('_')
    # Rejoin with -
    alphafold = '-'.join(alphafold)
    # Iterate through amyloid_dir that contains {amyloid}
    for amyloid_file in glob.glob(os.path.join(amyloid_dir,f'*{amyloid}*')):
        # Iterate through alphafold_dir that contains {alphafold}
        for alphafold_file in glob.glob(os.path.join(alphafold_dir,f'*{alphafold}*')):
            # Align the structures alphafold_file,amyloid_file,output_dir
            alignment(alphafold_file, amyloid_file, output_dir)
            
