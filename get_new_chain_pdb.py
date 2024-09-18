# get_new_chain_pdb.py
import os
import argparse

# Create an argument parser
parser = argparse.ArgumentParser(description='Extract a specific chain from a PDB file and rename it to chain A, while removing the original chain A.')

# Add optional arguments with '--' syntax
parser.add_argument('--pdb', required=True, help='Path to the input PDB file')
parser.add_argument('--chain', required=True, help='Chain identifier to extract (e.g., B)')

# Parse arguments
args = parser.parse_args()

# Get the base name of the file without extension
basename = os.path.splitext(os.path.basename(args.pdb))[0]

# Prepare the output file name
output_file = os.path.join(os.path.dirname(args.pdb), f"{basename}-chain-{args.chain}.pdb")

# Open the input PDB file for reading
with open(args.pdb, 'r') as infile:
    # Open the output PDB file for writing
    with open(output_file, 'w') as outfile:
        # Iterate through each line in the PDB file
        for line in infile:
            # Check if the line is an ATOM or HETATM record
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract the chain identifier (column 22)
                chain = line[21]
                # If the chain matches the specified chain ID, write it to the output file
                if chain == args.chain:
                    # Replace the chain ID with 'A' and write the modified line
                    modified_line = line[:21] + 'A' + line[22:]
                    outfile.write(modified_line)
                # Skip writing if the chain is the old chain A
                elif chain == 'A':
                    continue
            else:
                # Write non-ATOM/HETATM lines (e.g., headers, footers) to the output
                outfile.write(line)

print(f"New PDB file written to: {output_file}")
