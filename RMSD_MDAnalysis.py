#%%
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
import os
import pandas as pd
import numpy as np

# Define the paths to the FASTA folder and the reference PDB directory
fasta_dir = r"/mnt/c/Users/epicm/OneDrive/Desktop/Amyloid_Project/FASTA/alpha_synuclein_15aa_consecutive"
ref_dir = r"/mnt/c/Users/epicm/OneDrive/Desktop/Amyloid_Project/PDBs/alpha_synuclein_unique_amyloids"
segment_length = 15

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
    # Extract the ##-## part from the filename and return the first number
    numeric_part = filename.split('_')[2]
    start_num = int(numeric_part.split('-')[0])
    return start_num

#Create fasta_dict and pdb_dict
# Get a list of all rank_001 PDB files in the fasta directory
fasta_files = [f for f in os.listdir(fasta_dir) if "rank_001" in f and f.endswith(".pdb")]
# Sort the files based on the first number in the ##-## range
sorted_fasta_files = sorted(fasta_files, key=extract_range)
# Create a dictionary where the key is the first number in ##-## and the value is the file path
fasta_dict = {extract_range(f): os.path.join(fasta_dir, f) for f in sorted_fasta_files}
# Get a list of all pdb files in ref directory
ref_files = [f for f in os.listdir(ref_dir) if f.endswith(".pdb")]
# Create a dictionary where the key is the entry in the list starting at 1 and the value is the filepath
pdb_dict = {i: os.path.join(ref_dir, f) for i, f in enumerate(ref_files, start=1)}

# Initialize a dataframe to store RMSD values
# Determine the dimensions of the matrix
height = max(fasta_dict.keys()) + 1
length = len(pdb_dict)+1
#print height and length
print(f"Height: {height}, Length: {length}")
#initialize dataframe of height and length
df = pd.DataFrame(np.zeros((height, length)))
# Extract the base names without the .pdb extension
column_names = ["PDB ID"] + [os.path.splitext(os.path.basename(path))[0] for path in pdb_dict.values()]
# Set the column names
df.columns = column_names
# Set row names
df["PDB ID"] = [f"{i}:{i+14}" for i in range(1, height + 1)]

# Print the dictionary
#print(pdb_file_dict)
# Loop through each key in pdb_file_dict
    # Load ref universe pdb_file_dict[key]
    # Find from residue_id_dict first_resid and last_resid
        # Loop i from first_resid to last_resid - segment_length + 1
            # Load mobile universe fasta_dict[key], key = i
            # Select ref atoms as i:i+segment_length-1
            # Select mobile atoms as 1:15
            # Compute rmsd between mobile and ref
            # Store in df at row = i-1, column = key

# Loop through each key in pdb_file_dict
for key, ref_path in pdb_dict.items():

    # Set ref universe as ref_path
    ref_universe = mda.Universe(ref_path)

    # Find the first and last residue IDs for the reference structure
    first_resid, last_resid = find_chain_residues(ref_path)
    print(ref_path, first_resid, last_resid)

    # Loop from first_resid to last_resid-segment_length+1
    for i in range(first_resid, last_resid - segment_length + 1):

        # Set mobile universe to fasta_dict key = i
        mobile_path = fasta_dict[i]  # Assuming fasta_dict is defined and contains paths
        mobile_universe = mda.Universe(mobile_path)

        #select ref atoms from i to i+segment_length-1
        ref_atoms = ref_universe.select_atoms(f"chainID A and resid {i}:{i + segment_length - 1}")
        #select mobile atoms from 1:segment_length
        mobile_atoms = mobile_universe.select_atoms(f"chainID A and resid 1:{segment_length}")
        # Print ref_atoms vertically
        print("Reference Atoms:")
        for atom in ref_atoms:
            print(atom)

        # Print mobile_atoms vertically
        print("Mobile Atoms:")
        for atom in mobile_atoms:
            print(atom)
        
        # Reindex ref_atoms to start from 1


        # Compute RMSD between mobile and ref
        rmsd_value = rms.rmsd(mobile_atoms.positions, ref_atoms.positions, weights = "mass", center = True, superposition = True)
        print(f"RMSD for segment starting at residue {i}:{rmsd_value}")

        # Store RMSD in df at row = i-1, column = key
        df.loc[i, key] = rmsd_value
        
#%%
    

    

            









 # %%
# store as csv
output_csv = "rmsd_values.csv"
df.to_csv(output_csv, index=False)

#%%
# output the 5 lowest values in the RMSD 3 column
lowest_rmsd = df.nsmallest(5, "RMSD 3")

#%%
# histogram the RMSD 3 column
import matplotlib.pyplot as plt

plt.hist(df["RMSD 3"], bins=20, color='skyblue', edgecolor='black')
plt.xlabel("RMSD (nm)")
plt.ylabel("Frequency")