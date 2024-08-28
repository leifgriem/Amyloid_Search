# Initial_RMSD.py
# Calculate RMSD between AlphaFold predictions and reference PDB structures using MDAnalysis. Store results as a dataframe and optionally generate a heatmap.
#%%
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

# Define the paths to the FASTA folder and the reference PDB directory
fasta_dir = r"./PDBs/alpha_synuclein_15aa_consecutive_filtered"
ref_dir = r"./PDBs/alpha_synuclein_unique_amyloids"
output_dir = r"./PDBs/aligned_structures_initial"
# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)
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
    '''
    filename: str, the name of the pdb file
    '''
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
df = pd.DataFrame(np.full((height, length), np.inf))
# Extract the base names without the .pdb extension
column_names = ["PDB ID"] + [os.path.splitext(os.path.basename(path))[0] for path in pdb_dict.values()]
# Set the column names
df.columns = column_names
# Set row names
df["PDB ID"] = [f"{i}_{i+14}" for i in range(1, height + 1)]
# Create pdb_name_dict where the key is the same as in pdb_dict, and the value is the basename without the .pdb extension
pdb_name_dict = {key: os.path.splitext(os.path.basename(path))[0] for key, path in pdb_dict.items()}


# Loop through each key in pdb_file_dict
for key, ref_path in tqdm(pdb_dict.items()):

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

        # Define selection strings for mobile and reference universes
        mobile_selection = f"chainID A and resid {2*segment_length+21}:{3*segment_length+20} and name N CA C"
        ref_selection = f"chainID A and resid {i}:{i + segment_length - 1} and name N CA C"

        # Assign selections to ref_atoms and mobile_atoms
        ref_atoms = ref_universe.select_atoms(f"{ref_selection}")
        mobile_atoms = mobile_universe.select_atoms(f"{mobile_selection}")

        # Calculate RMSD
        rmsd_value = rms.rmsd(mobile_atoms.positions, ref_atoms.positions, center=True, superposition=True)
        # Convert from Angstroms to Nanometers
        rmsd_value = rmsd_value / 10
        df.loc[i, pdb_name_dict[key]] = rmsd_value

 # %%
import seaborn as sns
import matplotlib.pyplot as plt

# Create a mask for values less than or equal to 0
mask = df.iloc[:, 1:] <= 0  # Don't apply the mask to the first column with PDB IDs

# Set up the matplotlib figure
plt.figure(figsize=(10, 8))

# Sensitivity options
vmin = 0.1  # Set the minimum value for the colormap
vmax = 0.3  # Set the maximum value for the colormap

# Create a heatmap
sns.heatmap(
    df.iloc[:, 1:],  # Exclude the first column with PDB IDs
    mask=mask,  # Apply the mask
    cmap="viridis",  # Choose a colormap
    cbar_kws={"label": "RMSD (nm)"},  # Label for the colorbar
    annot=False,  # Don't show numbers in the cells
    vmin=vmin,  # Set minimum sensitivity
    vmax=vmax  # Set maximum sensitivity
)

# Set labels for the heatmap
plt.xlabel('Reference PDBs')
plt.ylabel('Alphafold Amyloid Predictions')

# Display the heatmap
plt.show()

# Save the heatmap as an image
plt.savefig("Initial_RMSD_Heatmap.png")

#%%
# Create a histogram of results
# %%
import matplotlib.pyplot as plt
import numpy as np

# Filter out np.inf values to prepare data for histogram
rmsd_values = df.iloc[:, 1:].to_numpy().flatten()  # Flatten the DataFrame to a 1D array
rmsd_values = rmsd_values[np.isfinite(rmsd_values)]  # Remove np.inf values

# Create a histogram of RMSD values
plt.figure(figsize=(8, 6))
plt.hist(rmsd_values, bins=30, color='blue', alpha=0.7, edgecolor='black')
plt.title('Histogram of RMSD Values')
plt.xlabel('RMSD (nm)')
plt.ylabel('Frequency')

# Display the histogram
plt.show()

# Save the histogram as an image
plt.savefig("Initial_RMSD_Histogram.png")



# %%
