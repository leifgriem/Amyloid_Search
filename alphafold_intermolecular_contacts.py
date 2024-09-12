# alphafold_intermolecular_contacts.py
# Generates a csv and histogram of intermolecular contacts for AlphaFold fragments
import argparse
import numpy as np
import pandas as pd
from glob import glob
import mdtraj as md
import os
from itertools import combinations
from matplotlib import pyplot as plt
from tqdm import tqdm
import re

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Analyze AlphaFold PDB files and generate CSV and histogram outputs.")

# Add arguments for directory paths, protein name, and segment length
parser.add_argument('--alphafold_directory', type=str, required=True, help='Directory containing AlphaFold PDB files.')
parser.add_argument('--csv_output_dir', type=str, required=True, help='Directory to save CSV output.')
parser.add_argument('--histogram_output_dir', type=str, required=True, help='Directory to save histogram images.')
parser.add_argument('--protein', type=str, required=True, help='Name of the protein being analyzed.')
parser.add_argument('--segment_length', type=int, required=True, help='Segment length.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
alphafold_directory = args.alphafold_directory
csv_output_dir = args.csv_output_dir
histogram_output_dir = args.histogram_output_dir
protein = args.protein
segment_length = args.segment_length

# Ensure output directories exist
os.makedirs(csv_output_dir, exist_ok=True)
os.makedirs(histogram_output_dir, exist_ok=True)

def get_pairdistances(peptide_path):
    AF2_pdb = md.load(peptide_path)
    topology = AF2_pdb.topology

    alpha_carbons = ([atom.index for atom in topology.atoms if atom.name == 'CA'])

    atom_pairs = list(combinations(alpha_carbons, 2))
    pairwise_distances = md.geometry.compute_distances(AF2_pdb, atom_pairs)[0]

    num_residues = AF2_pdb.n_residues  # Number of residues in the protein
    atom_to_residue = {atom.index: atom.residue.index for atom in topology.atoms if atom.name == 'CA'}

    # Initialize an empty 2D matrix for the distances
    distance_matrix = np.zeros((num_residues, num_residues))

    # Fill the distance matrix. Since the distances are symmetric, we mirror the values across the diagonal
    for distance, (atom_index_i, atom_index_j) in zip(pairwise_distances, atom_pairs):
        residue_index_i = atom_to_residue[atom_index_i]
        residue_index_j = atom_to_residue[atom_index_j]

        # Populate the matrix, adjusting indices by -1 if necessary
        # Adjust the indices based on how your residues are indexed (0-based or 1-based)
        distance_matrix[residue_index_i][residue_index_j] = distance
        distance_matrix[residue_index_j][residue_index_i] = distance  # Mirror the distance for symmetry
    assert distance_matrix.shape == (segment_length * 5, segment_length * 5), f'Expected 75x75 array and got {np.size(distance_matrix)}'

    # Convert to pandas df
    chains = ['A', 'B', 'C', 'D', 'E']
    resids = np.arange(1, segment_length + 1)
    indices = [f'{chain}{resid}' for chain in chains for resid in resids]
    distance_df = pd.DataFrame(distance_matrix, columns=indices, index=indices)
    return distance_df

def get_chain_info(peptide_path):
    AF2_pdb = md.load(peptide_path)
    topology = AF2_pdb.topology

    if peptide_path.split("/")[0] in ['multimer', 'multimer_5rec']:
        chain_info = []
        for chain in topology.chains:
            chain_residues = [residue.index for residue in chain.residues]
            chain_start = min(chain_residues)
            chain_end = max(chain_residues)
            chain_info.append((chain_start, chain_end))
    # Hard coded chain info for 10U since every residue is chain A
    else:
        chain_info = [(0, segment_length - 1), (segment_length, 2 * segment_length - 1), 
                      (2 * segment_length, 3 * segment_length - 1), (3 * segment_length, 4 * segment_length - 1), 
                      (4 * segment_length, 5 * segment_length - 1)]

    return chain_info

def count_interchain_contacts(distance_matrix, peptide_path):
    contact_map = (distance_matrix < 0.8).values  # Assuming <0.8 nm as contact threshold
    interchain_contacts = 0
    chain_info = get_chain_info(peptide_path)
    for i, (start_i, end_i) in enumerate(chain_info):
        for j, (start_j, end_j) in enumerate(chain_info):
            if i < j:  # Avoid double counting and self-counting
                # Slice the contact map to only consider contacts between chains i and j
                interchain_section = contact_map[start_i:end_i + 1, start_j:end_j + 1]
                interchain_contacts += np.sum(interchain_section)

    return interchain_contacts

# List all the pdb files in the directory
alphafold_pdbs = glob(f'{alphafold_directory}/*.pdb')

# Sort the pdbs by the fragment number
def extract_first_number(filename):
    # Split the filename to extract the range part
    parts = filename.split('_')
    for part in parts:
        if re.match(r'^\d+-\d+$', part):
            return int(part.split('-')[0])

alphafold_pdbs = sorted(alphafold_pdbs, key=extract_first_number)

# Get just the pdb number
alphafold_names = [pdb.split('/')[-1].split('_unrelaxed')[0] for pdb in alphafold_pdbs]
alphafold_names = [name.split('_')[-1] for name in alphafold_names]

# Calculate interchain contacts
interchain_contacts = []
for pdb, name in tqdm(zip(alphafold_pdbs, alphafold_names), desc="Processing PDB files", total=len(alphafold_pdbs)):
    distance_df = get_pairdistances(pdb)
    interchain_contacts.append(count_interchain_contacts(distance_df, pdb))

# normalize the interchain contacts
normed_interchain_contacts = [contacts / (segment_length **2 ) for contacts in interchain_contacts]

# Output a CSV with the interchain contacts
interchain_contacts_df = pd.DataFrame(normed_interchain_contacts, columns=['Interchain Contacts'])
interchain_contacts_df.index = alphafold_names
csv_output_path = os.path.join(csv_output_dir, f'normalized_interchain_contacts_{protein}_{segment_length}aa.csv')
interchain_contacts_df.to_csv(csv_output_path, index=True)
print("CSV output successfully saved.")

# Plot a histogram of the interchain contacts
plt.hist(normed_interchain_contacts, bins=15)
plt.xlabel('Number of Interchain Contacts')
plt.ylabel('Number of Fragments')
plt.title(f"Normalized AlphaFold Fragments' Interchain Contacts: {protein}_{segment_length}aa")

# Save histogram as an image
histogram_output_path = os.path.join(histogram_output_dir, f'normalized_interchain_contacts_{protein}_{segment_length}aa.png')
plt.savefig(histogram_output_path)
# plt.show()