#run_analysis.py
#%%
import subprocess
import os

# For 13aa
# alphafold_raw_output = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/alpha_synuclein_fragments_13aa"  # Directory containing raw AlphaFold output
# alphafold_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/alpha-synuclein_13aa_consecutive_filtered"
# ref_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/alpha-synuclein_unique_amyloids"
# csv_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/CSVs"
# img_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/IMGs"
# aligned_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/aligned_structures_13aa"  # Output directory for Alignment_MDAnalysis
# aligned_reindexed_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/aligned_reindexed_structures_13aa"  # Output directory for fix_resid_chainid
# protein = "alpha-synuclein"
# segment_length = 13
num = 19
# Define directories and common arguments
alphafold_raw_output = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/alpha_synuclein_fragments_{num}aa"  # Directory containing raw AlphaFold output
alphafold_dir = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/alpha-synuclein_{num}aa_consecutive_filtered"
ref_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/alpha-synuclein_unique_amyloids"
csv_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/CSVs"
img_dir = "/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/IMGs"
aligned_dir = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/aligned_structures_{num}aa"  # Output directory for Alignment_MDAnalysis
aligned_reindexed_dir = f"/mnt/c/Users/epicm/OneDrive/Desktop/amyloid_search/alpha-synuclein/PDBs/aligned_reindexed_structures_{num}aa"  # Output directory for fix_resid_chainid
protein = "alpha-synuclein"
segment_length = num

def run_script(script_name, args):
    try:
        print(f"Running {script_name}...")
        result = subprocess.run(['python', script_name] + args, check=True, text=True, capture_output=False)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e.stderr}")
        exit(1)
#%%
# Run get_rank_001.py to copy rank_001 PDB files
run_script("get_rank_001.py", [
    "--alphafold_raw_output", alphafold_raw_output,
    "--output_dir", alphafold_dir
])

#%%
# Run Initial RMSD calculation
run_script("Initial_RMSD.py", [
    "--alphafold_dir", alphafold_dir,
    "--ref_dir", ref_dir,
    "--csv_dir", csv_dir,
    "--img_dir", img_dir,
    "--protein", protein,
    "--segment_length", str(segment_length)
])
#%%
# Run AlphaFold intermolecular contacts calculation
run_script("alphafold_intermolecular_contacts.py", [
    "--alphafold_directory", alphafold_dir,
    "--csv_output_dir", csv_dir,
    "--histogram_output_dir", img_dir,
    "--protein", protein,
    "--segment_length", str(segment_length)
])
#%%
# Run Alignment MDAnalysis
run_script("Alignment_MDAnalysis.py", [
    "--interchain_contacts", os.path.join(csv_dir, f"normalized_interchain_contacts_{protein}_{segment_length}aa.csv"),
    "--RMSDs", os.path.join(csv_dir, f"initial_RMSD_{protein}_{segment_length}aa.csv"),
    "--alphafold_dir", alphafold_dir,
    "--amyloid_dir", ref_dir,
    "--output_dir", aligned_dir,  # Output directory for aligned structures
    "--protein", protein,
    "--segment_length", str(segment_length)
])
#%%
# Run fix_resid_chainid (takes aligned_dir as input and outputs to aligned_reindexed_dir)
run_script("fix_resid_chainid.py", [
    "--input_dir", aligned_dir,  # Input directory from Alignment_MDAnalysis output
    "--output_dir", aligned_reindexed_dir,  # Output directory for reindexed structures
    "--segment_length", str(segment_length),
    "--protein", protein
])
#%%
# Run RMSD after filtering (takes aligned_reindexed_dir as input)
run_script("rmsd_after_filtering.py", [
    "--alphafold_dir", aligned_reindexed_dir,  # Input directory from fix_resid_chainid output
    "--amyloid_dir", ref_dir,
    "--csv_dir", csv_dir,
    "--img_dir", img_dir,
    "--protein", protein,
    "--segment_length", str(segment_length)
])

print("All tasks completed.")

# %%
