# Chimera Processing
#%% # Fix sequencing issue in my FASTA file
import re

def flip_sequence(sequence):
    """Flip the amino acid sequence before PGGG and after PGGG, ignoring U's"""
    if "PGGG" in sequence:
        parts = sequence.split("PGGG")
        # Flip the first part normally
        before_pggg = parts[0][::-1]
        # For the second part, exclude U's from the flip
        after_pggg = parts[1].rstrip('U')[::-1] + 'U' * parts[1].count('U')
        return before_pggg + "PGGG" + after_pggg
    else:
        return sequence

def process_fasta(fasta_file, output_file):
    with open(fasta_file, "r") as infile, open(output_file, "w") as outfile:
        buffer = []
        header = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                # Process previous buffer if there was one
                if header:
                    write_entry(header, buffer, outfile)
                # Start new entry
                header = line
                buffer = []
            else:
                # Collect sequence lines (up to 5)
                buffer.append(line)
                if len(buffer) == 5:
                    # Automatically process once 5 lines are collected
                    write_entry(header, buffer, outfile)
                    header = ""
                    buffer = []

        # Process any leftover entry in the buffer
        if header:
            write_entry(header, buffer, outfile)

def write_entry(header, sequence_lines, outfile):
    """Process and write the entry to the output file"""
    # Extract numbering from the header
    match = re.match(r'>(\w+)_([0-9]+)-([0-9]+)_PGGG_([0-9]+)-([0-9]+)', header)
    if match:
        pdb_name = match.group(1)
        num1, num2, num3, num4 = int(match.group(2)), int(match.group(3)), int(match.group(4)), int(match.group(5))

        if num1 > num2:
            # Flip the numbers and sequences
            new_header = f">{pdb_name}_{num2}-{num1}_PGGG_{num4}-{num3}\n"
            flipped_sequences = [flip_sequence(seq) for seq in sequence_lines]
            outfile.write(new_header)
            for seq in flipped_sequences:
                outfile.write(seq + '\n')
        else:
            # Write the original header and sequences
            outfile.write(header + '\n')
            for seq in sequence_lines:
                outfile.write(seq + '\n')

# Example usage
input_fasta = "./FASTA/alpha_synuclein_chimeras.fasta"
output_fasta = "./FASTA/alpha_synuclein_chimeras_renumbered.fasta"
process_fasta(input_fasta, output_fasta)

#%% # Get rank 001 files from the raw output and put them in new directory
import os
import shutil

# Define the source and destination directories
source_dir = './PDBs/alpha_synuclein_chimeras_refined'
dest_dir = './PDBs/alpha_synuclein_chimeras_refined_rank_001'

# Ensure the destination directory exists, if not, create it
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

# Iterate over all files in the source directory
for filename in os.listdir(source_dir):
    # Check if the file contains 'rank_001' and ends with '.pdb'
    if 'rank_001' in filename and filename.endswith('.pdb'):
        # Construct the full file paths
        src_file = os.path.join(source_dir, filename)
        dest_file = os.path.join(dest_dir, filename)
        
        # Copy the file to the destination directory
        shutil.copy(src_file, dest_file)
        print(f"Copied: {filename}")


#%% # Reindex chimera PDB files
import os
from glob import glob

# Specify the input directory containing the PDB files and output directory
input_dir = './PDBs/alpha_synuclein_chimeras_refined_rank_001'  # Replace with your input directory path containing PDB files
output_dir = './PDBs/chimera_pdbs_refined_renumbered'  # Output directory for reindexed PDB files

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def fix_resid_chainid_chimera(pdb_path, output_dir):
    # Extract just the base name of the PDB file
    pdb_base_name = os.path.basename(pdb_path)
    
    # Read the input PDB file
    with open(pdb_path, 'r') as file:
        lines = file.readlines()

    # Count the total number of atoms
    total_atoms = sum(1 for line in lines if line.startswith(("ATOM", "HETATM")))

    # Calculate the number of atoms per chain
    atoms_per_chain = total_atoms // 5
    remainder_atoms = total_atoms % 5  # To handle any remainder atoms

    # Initialize variables for chain and residue numbering
    chain_assignment = ['A', 'B', 'C', 'D', 'E']
    current_chain_index = 0
    current_atom_count = 0
    current_resid = 0  # Start at 0 to handle initial increment correctly

    # Track the last processed residue ID to detect changes
    last_resid = None

    output_lines = []

    # Iterate over each line to assign atoms to chains and renumber residues
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            current_atom_count += 1
            resid = int(line[22:26].strip())

            # Check if we need to move to the next chain
            if (current_atom_count > atoms_per_chain and current_chain_index < len(chain_assignment) - 1) or (current_atom_count == atoms_per_chain + 1 and remainder_atoms > 0):
                current_chain_index += 1
                current_atom_count = 1  # Reset atom count for new chain
                current_resid = 0  # Reset residue number for new chain
                remainder_atoms -= 1  # Assign extra atoms to the first few chains
                last_resid = None  # Reset last_resid for the new chain

            current_chain = chain_assignment[current_chain_index]

            # Increment the residue number only if the residue ID has changed
            if last_resid is None or resid != last_resid:
                current_resid += 1
                last_resid = resid

            # Update line with new chain ID and residue ID
            line = line[:21] + current_chain + line[22:]
            line = line[:22] + f"{current_resid:4d}" + line[26:]

        output_lines.append(line)

    # Define the output path for the reindexed file
    output_path = os.path.join(output_dir, pdb_base_name)

    # Write the modified lines to the new PDB file
    with open(output_path, 'w') as file:
        file.writelines(output_lines)

# Loop through all PDB files in the input directory and apply chain renumbering
pdb_files = glob(os.path.join(input_dir, '*.pdb'))
for pdb_file in pdb_files:
    fix_resid_chainid_chimera(pdb_file, output_dir)




#%% # Find normalized interchain contacts and copy to chimera_pdbs_filtered
import os
import shutil
import numpy as np
import mdtraj as md
from itertools import combinations
from matplotlib import pyplot as plt
from glob import glob

# Specify the directory containing your PDB files
pdb_directory = './PDBs/chimera_pdbs_refined_renumbered'
filtered_pdb_directory = './PDBs/chimera_pdbs_filtered/'

# Create the new directory if it doesn't exist
os.makedirs(filtered_pdb_directory, exist_ok=True)

def get_pairdistances(peptide_path):
    AF2_pdb = md.load(peptide_path)
    topology = AF2_pdb.topology

    # Get alpha carbon (CA) atoms
    alpha_carbons = [atom.index for atom in topology.atoms if atom.name == 'CA']
    
    # Get all pairwise combinations of alpha carbon atoms
    atom_pairs = list(combinations(alpha_carbons, 2))
    
    # Compute pairwise distances between alpha carbons
    pairwise_distances = md.geometry.compute_distances(AF2_pdb, atom_pairs)[0]

    num_residues = AF2_pdb.n_residues
    atom_to_residue = {atom.index: atom.residue.index for atom in topology.atoms if atom.name == 'CA'}

    # Initialize distance matrix
    distance_matrix = np.zeros((num_residues, num_residues))

    # Fill the distance matrix with computed distances
    for distance, (atom_index_i, atom_index_j) in zip(pairwise_distances, atom_pairs):
        residue_index_i = atom_to_residue[atom_index_i]
        residue_index_j = atom_to_residue[atom_index_j]

        distance_matrix[residue_index_i][residue_index_j] = distance
        distance_matrix[residue_index_j][residue_index_i] = distance

    return distance_matrix

def count_interchain_contacts(distance_matrix, peptide_path):
    # Create contact map by checking distances < 1.5 nm (or any threshold you prefer)
    contact_map = (distance_matrix < 1.5).astype(int)
    interchain_contacts = 0
    chain_info = get_chain_info(peptide_path)
    
    # Create a matrix to visualize interchain contacts
    contact_map_matrix = np.zeros_like(distance_matrix)

    # Iterate through the chain pairs
    for i, (start_i, end_i) in enumerate(chain_info):
        for j, (start_j, end_j) in enumerate(chain_info):
            if i != j:  # Handle both interchain contacts (i < j and i > j)
                # Get interchain contact section
                interchain_section = contact_map[start_i:end_i + 1, start_j:end_j + 1]
                interchain_contacts += np.sum(interchain_section)

                # Populate both upper and lower triangular sections
                contact_map_matrix[start_i:end_i + 1, start_j:end_j + 1] = interchain_section
                contact_map_matrix[start_j:end_j + 1, start_i:end_i + 1] = interchain_section  # Mirror to lower triangle

    # Calculate the segment length for normalization (assuming both chains have the same length)
    segment_length = min(end_i - start_i + 1, end_j - start_j + 1)  # Use the minimum segment length
    normalized_interchain_contacts = interchain_contacts / (segment_length ** 2)

    return interchain_contacts, contact_map_matrix, normalized_interchain_contacts

def get_chain_info(peptide_path):
    AF2_pdb = md.load(peptide_path)
    topology = AF2_pdb.topology

    # Extract chain start and end residues
    chain_info = []
    for chain in topology.chains:
        chain_residues = [residue.index for residue in chain.residues]
        chain_start = min(chain_residues)
        chain_end = max(chain_residues)
        chain_info.append((chain_start, chain_end))

    return chain_info

# List all PDB files in the directory
pdb_files = glob(os.path.join(pdb_directory, '*.pdb'))

# Store normalized interchain contacts for each PDB file
normalized_contacts_list = []

# Loop through all PDB files and calculate interchain contacts
for pdb_file in pdb_files:
    distance_matrix = get_pairdistances(pdb_file)
    _, _, normalized_interchain_contacts = count_interchain_contacts(distance_matrix, pdb_file)

    # Add the normalized interchain contacts to the list
    normalized_contacts_list.append((pdb_file, normalized_interchain_contacts))

# Calculate the average value of normalized interchain contacts
average_contacts = np.mean([norm_contacts for _, norm_contacts in normalized_contacts_list])

# Filter PDB files with normalized interchain contacts above the average and copy to the new directory
for pdb_file, norm_contacts in normalized_contacts_list:
    if norm_contacts > average_contacts:
        # Copy the file to the filtered directory
        shutil.copy(pdb_file, os.path.join(filtered_pdb_directory, os.path.basename(pdb_file)))

# Plot a histogram of the normalized interchain contacts for all PDB files
plt.figure(figsize=(8, 6))
plt.hist([norm_contacts for _, norm_contacts in normalized_contacts_list], bins=15, color='blue', alpha=0.7, edgecolor='black')
plt.axvline(average_contacts, color='red', linestyle='dashed', linewidth=1, label=f'Average = {average_contacts:.2f}')
plt.xlabel('Normalized Interchain Contacts')
plt.ylabel('Frequency')
plt.title(f"Normalized Interchain Contacts for All PDB Files in {pdb_directory}")
plt.legend()
plt.show()

# Print the average normalized interchain contacts
print(f"Average normalized interchain contacts: {average_contacts:.3f}")



#%% # Align chimera PDBs to amyloid PDBs
import os
import re
import MDAnalysis as mda
from MDAnalysis.analysis import align

# Directories for input PDB files and output aligned structures
filtered_dir = './PDBs/chimera_pdbs_refined_renumbered'
amyloid_dir = './PDBs/alpha-synuclein_unique_amyloids'
output_dir = './PDBs/aligned_refined_chimera_onto_amyloids'

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

def extract_residue_range_and_code(filename):
    """
    Extracts num1, num2, num3, num4, and PDB_name from the file path by splitting on underscores.
    """
    # Get the filename without the directory and extension
    filename = os.path.basename(filename).replace(".pdb", "")
    
    # Split the filename by underscores
    parts = filename.split("_")
    
    # Extract values based on specified indexes
    PDB_name = parts[0]
    num1, num2 = map(int, parts[1].split("-"))
    num3, num4 = map(int, parts[3].split("-"))
    
    return PDB_name, num1, num2, num3, num4

def align_and_save_structure(alphafold_universe, amyloid_universe, filtered_resid_selection, amyloid_resid_selection, output_dir, filtered_pdb, amyloid_name):
    # Perform MDAnalysis align.alignto command using the backbone selection
    align.alignto(alphafold_universe, amyloid_universe, select=(filtered_resid_selection, amyloid_resid_selection), weights="mass")

    # After alignment, we want to write the full aligned structure (not just the selected atoms)
    all_atoms = alphafold_universe.atoms

    # Set the filename for the aligned structure
    aligned_filename = f"{os.path.splitext(filtered_pdb)[0]}_{amyloid_name}_aligned.pdb"
    
    # Save the entire aligned structure to the output directory
    aligned_path = os.path.join(output_dir, aligned_filename)
    
    # Write the entire aligned structure to the output directory
    with mda.Writer(aligned_path, all_atoms.n_atoms) as W:
        W.write(all_atoms)

    print(f"Aligned structure saved as: {aligned_path}")

# Loop through all PDB files in filtered_dir
for filtered_pdb in os.listdir(filtered_dir):
    if filtered_pdb.endswith('.pdb'):
        print(f"Processing filtered PDB: {filtered_pdb}")  # Print statement for the filtered PDB being processed

        # Load the filtered PDB structure
        filtered_pdb_path = os.path.join(filtered_dir, filtered_pdb)
        alphafold_universe = mda.Universe(filtered_pdb_path)

        # Extract num1, num2, num3, and num4 from the filename
        PDB_name, num1, num2, num3, num4 = extract_residue_range_and_code(filtered_pdb)

        # Define the corresponding amyloid PDB filename
        amyloid_pdb = f"{PDB_name}.pdb"  # Use the string in position 0 as the amyloid PDB name
        amyloid_pdb_path = os.path.join(amyloid_dir, amyloid_pdb)

        # Check if the corresponding amyloid PDB exists
        if not os.path.exists(amyloid_pdb_path):
            print(f"Amyloid PDB {amyloid_pdb} not found, skipping.")
            continue

        print(f"Processing amyloid PDB: {amyloid_pdb}")  # Print statement for the amyloid PDB being processed

        # Load the amyloid PDB structure
        amyloid_universe = mda.Universe(amyloid_pdb_path)

        # Calculate num_diff and num_diff2
        num_diff = abs(num2 - num1) + 1
        num_diff2 = abs(num4 - num3) + 1

        # Define residue selection for the filtered PDB (chain C)
        # Select residues 1:num_diff and num_diff+5:num_diff+5+num_diff2
        filtered_resid_selection = f"chainID C and (resid 1:{num_diff} or resid {num_diff+5}:{num_diff+4+num_diff2}) and name N CA C"

        # Print the filtered_resid_selection to verify correctness
        print(f"Filtered residue selection: {filtered_resid_selection}")

        # Special case for PDB 6sst_85-91_PGGG_15-20_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb
        if filtered_pdb == '6sst_85-91_PGGG_15-20_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb':
            amyloid_resid_selection = "(chainID A and resid 85:91) or (chainID B and resid 15:20) and name N CA C"
        # Handle interchain chimera case: if num1 = num3 and num2 = num4, select from chain A and chain B
        elif num1 == num3 and num2 == num4:
            amyloid_resid_selection = f"(chainID A and resid {num1}:{num2}) or (chainID B and resid {num3}:{num4}) and name N CA C"
        else:
            amyloid_resid_selection = f"chainID A and (resid {num1}:{num2} or resid {num3}:{num4}) and name N CA C"

        # Print the amyloid_resid_selection to verify correctness
        print(f"Amyloid residue selection: {amyloid_resid_selection}")

        # Align and save the entire structure
        amyloid_name = os.path.splitext(amyloid_pdb)[0]  # Use the amyloid PDB name
        align_and_save_structure(alphafold_universe, amyloid_universe, filtered_resid_selection, amyloid_resid_selection, output_dir, filtered_pdb, amyloid_name)

# %% # Rearrange chains for visualization
import os

# Directory with the aligned PDB files
pdb_dir = './PDBs/aligned_refined_chimera_onto_amyloids'
output_dir = './PDBs/chimeras_refined_visual'

# Make sure the output directory exists
os.makedirs(output_dir, exist_ok=True)
# Make the output directory if it doesnt exist


# Mapping chain order: B -> B, C -> A, D -> C
chain_mapping = {'B': 'B', 'C': 'A', 'D': 'C'}

# Helper function to extract the simplified filename format
def extract_new_filename(pdb_file):
    parts = pdb_file.split('_')
    # The new format: PDBname_num1-num2_PGGG_num3-num4.pdb
    return f"{parts[0]}_{parts[1]}_{parts[2]}_{parts[3]}.pdb"

# Loop through all PDB files in the input directory
for pdb_file in os.listdir(pdb_dir):
    if pdb_file.endswith('.pdb'):
        print(f"Processing {pdb_file}...")

        # Full path to the input PDB file
        pdb_path = os.path.join(pdb_dir, pdb_file)

        # Create the new output filename based on the simplified format
        new_filename = extract_new_filename(pdb_file)

        # Full path for the output file
        output_file = os.path.join(output_dir, new_filename)

        # Process the PDB file
        with open(pdb_path, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                # Only modify lines starting with ATOM or HETATM
                if line.startswith(('ATOM', 'HETATM')):
                    chain_id = line[21]  # Chain identifier is at position 22 (index 21)

                    # Skip chain A or E if we are deleting them
                    if chain_id in ['A', 'E']:
                        continue

                    # Reindex chains according to the chain_mapping dictionary
                    if chain_id in chain_mapping:
                        new_chain_id = chain_mapping[chain_id]
                        # Replace the chain ID with the new chain ID in the correct column
                        line = line[:21] + new_chain_id + line[22:]

                # Write the modified or unmodified line to the output PDB
                outfile.write(line)

        print(f"Rearranged PDB saved as: {output_file}")

print("All PDB files have been processed.")

# %%
