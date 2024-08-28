#alphafold_intermolecular_contacts.py
#Generates a csv and histogram of intermolecular contacts for AlphaFold fragments
#%%
import numpy as np
import pandas as pd
from glob import glob
import mdtraj as md
import os
from itertools import combinations
from matplotlib import pyplot as plt
from tqdm import tqdm

protein = 'alpha_synuclein'

def get_pairdistances(peptide_path):
    AF2_pdb = md.load(peptide_path)
    topology = AF2_pdb.topology

    alpha_carbons = ([atom.index for atom in topology.atoms if atom.name == 'CA'])

    atom_pairs = list(combinations(alpha_carbons,2))
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
        # print(distance_matrix)
        # print(np.shape(distance_matrix))
    assert distance_matrix.shape == (75, 75), f'Expected 75x75 array and got {np.size(distance_matrix)}'

    # convert to pandas df
    chains = ['A', 'B', 'C', 'D', 'E']
    resids = np.arange(1,16)
    indices = [f'{chain}{resid}' for chain in chains for resid in resids]
    distance_df = pd.DataFrame(distance_matrix, columns=indices, index=indices)
    return distance_df

def get_chain_info(peptide_path):
    AF2_pdb = md.load(peptide_path)
    topology = AF2_pdb.topology
    
    if peptide_path.split("/")[0] in ['multimer','multimer_5rec']:
        chain_info = []
        for chain in topology.chains:
            chain_residues = [residue.index for residue in chain.residues]
            chain_start = min(chain_residues)
            chain_end = max(chain_residues)
            chain_info.append((chain_start, chain_end))
    # hard coded chain info for 10U since every residue is chain A
    else: chain_info = [(0, 14), (15, 29), (30, 44), (45, 59), (60, 74)]

    return chain_info

def count_interchain_contacts(distance_matrix, peptide_path):
    contact_map = (distance_matrix < 0.8).values  # Assuming <0.8 nm as contact threshold
    interchain_contacts = 0
    chain_info = get_chain_info(peptide_path)
    print(chain_info)
    for i, (start_i, end_i) in enumerate(chain_info):
        for j, (start_j, end_j) in enumerate(chain_info):
            if i < j:  # Avoid double counting and self-counting
                # Slice the contact map to only consider contacts between chains i and j
                interchain_section = contact_map[start_i:end_i+1, start_j:end_j+1]
                # print(interchain_section)
                interchain_contacts += np.sum(interchain_section)
    
    return interchain_contacts

def plot_pairdistances(distance_matrix, peptide_path, peptide_name):
    AF2_pdb = md.load(peptide_path)
    topology = AF2_pdb.topology
    chain_info = get_chain_info(peptide_path)
    chain_starts = [start for start, _ in chain_info]

    plt.figure(figsize=(10, 8))
    plt.imshow(distance_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Distance (nm)')
    
    # Draw lines to demarcate different chains
    for start in chain_starts:
        # print(start)
        plt.axvline(x=start-0.5, color='black')
        plt.axhline(y=start-0.5, color='black')

    plt.title(f'Pairwise Distance Heatmap: {peptide_name}')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    # plt.show()

def get_contactmap(distance_matrix, peptide_path, peptide_name):

    # output mask of contacts (<0.6nm)
    contact_map = distance_matrix < 0.6
    chain_info = get_chain_info(peptide_path)
    chain_starts = [start for start, _ in chain_info]

    plt.figure(figsize=(10, 8))
    plt.imshow(contact_map, cmap='viridis', interpolation='nearest')
    
    # Draw lines to demarcate different chains
    for start in chain_starts:
        plt.axvline(x=start-0.5, color='black')
        plt.axhline(y=start-0.5, color='black')

    plt.title(f'Contact Map: {peptide_name}')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    # plt.show()

#%%

alphafold_directory = f'PDBs/{protein}_15aa_consecutive_filtered'
# list all the pdb files in the directory
alphafold_pdbs = glob(f'{alphafold_directory}/*.pdb')
# sort the pdbs by the fragment number
def extract_first_number(filename):
    # Split the filename to extract the range part
    parts = filename.split('_')
    for part in parts:
        if '-' in part:
            return int(part.split('-')[0])
alphafold_pdbs = sorted(alphafold_pdbs, key=extract_first_number)

# get just the pdb number
# split the path by / and get the last part, then split by _unrelaxed and get the first part
alphafold_names = [pdb.split('/')[-1].split('_unrelaxed')[0] for pdb in alphafold_pdbs]
# Split by _ and get the last part
alphafold_names = [name.split('_')[-1] for name in alphafold_names]

#%%
interchain_contacts = []
for pdb,name in tqdm(zip(alphafold_pdbs, alphafold_names)):
    distance_df = get_pairdistances(pdb)
    plot_pairdistances(distance_df, pdb, name)
    get_contactmap(distance_df, pdb, name)
    interchain_contacts.append(count_interchain_contacts(distance_df, pdb))

#%%
# output a csv with the interchain contacts
interchain_contacts_df = pd.DataFrame(interchain_contacts, columns=['Interchain Contacts'])
interchain_contacts_df.index = alphafold_names
interchain_contacts_df.to_csv(f'interchain_contacts_{protein}.csv', index=True)
print("outputted csv successfully")
print(interchain_contacts_df)

#%%
# plot a histogram of the interchain contacts
plt.hist(interchain_contacts, bins=15)
plt.xlabel('Number of Interchain Contacts')
plt.ylabel('Number of Fragments')
plt.title(f"AlphaFold Fragments' Interchain Contacts: {protein}")
# %%
