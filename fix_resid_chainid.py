import argparse
import os
from tqdm import tqdm

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Reindex PDB files based on segment length and save them to a specified directory.")

# Add arguments for input directory, output directory, segment length, and protein name
parser.add_argument('--input_dir', type=str, required=True, help='Directory containing the input PDB files.')
parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the reindexed PDB files.')
parser.add_argument('--segment_length', type=int, required=True, help='Segment length.')
parser.add_argument('--protein', type=str, required=True, help='Name of the protein being analyzed.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
input_dir = args.input_dir
output_dir = args.output_dir
segment_length = args.segment_length
protein = args.protein

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def fix_resid_chainid(pdb, output_dir):

    # Extract just the base name of the PDB file
    pdb_base_name = os.path.basename(pdb)
    
    # Extract range from the input PDB filename
    pdb_name = pdb_base_name.split('_')
    pdb_resid = pdb_name[1]
    first_resid, last_resid = map(int, pdb_resid.split('-'))

    # This dictionary maps original residue ranges to new residue IDs and chains
    residue_mapping = {
        range(segment_length + 11, 2 * segment_length + 11): (range(first_resid, last_resid + 1), 'B'),  # First segment to chain B
        range(2 * segment_length + 21, 3 * segment_length + 21): (range(first_resid, last_resid + 1), 'A'),  # Middle segment to chain A
        range(3 * segment_length + 31, 4 * segment_length + 31): (range(first_resid, last_resid + 1), 'C'),  # Last segment to chain C
    }

    # Read the input PDB file
    with open(pdb, 'r') as file:
        lines = file.readlines()

    output_lines = []

    # Iterate through each line
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Extract the current residue ID and chain ID
            resid = int(line[22:26].strip())
            chain = line[21].strip()

            # Check if the residue ID falls within any of the defined ranges
            for resid_range, (new_resids, new_chain) in residue_mapping.items():
                if resid in resid_range:
                    # Calculate the new residue ID based on the index
                    new_resid = new_resids[resid - resid_range.start]

                    # Replace the chain ID and residue ID in the current line
                    line = line[:21] + new_chain + line[22:]
                    line = line[:22] + f"{new_resid:4d}" + line[26:]

                    break  # No need to check further, as the line has been updated

        # Add the modified or unmodified line to the output list
        output_lines.append(line)

    # Define the output path for the reindexed file
    output_path = os.path.join(output_dir, pdb_base_name)

    # Write the modified lines to the new PDB file
    with open(output_path, 'w') as file:
        file.writelines(output_lines)

# Loop through all PDB files in the input directory and run the function
pdb_files = [file for file in os.listdir(input_dir) if file.endswith('.pdb')]
for file in tqdm(pdb_files, desc="Processing PDB files"):
    input_path = os.path.join(input_dir, file)
    fix_resid_chainid(input_path, output_dir)

