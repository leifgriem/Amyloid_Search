import argparse
import os

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Process amino acid sequence, fragment length, protein name, and output directory.")

# Add arguments for sequence, protein name, and output directory
parser.add_argument('--sequence', type=str, required=True, help='Amino acid sequence to be processed.')
parser.add_argument('--protein_name', type=str, required=True, help='Name of the protein.')
parser.add_argument('--output_dir', type=str, required=True, help='Directory where output files will be saved.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
sequence = args.sequence
protein_name = args.protein_name
output_dir = args.output_dir

# Fragment lengths to be processed (only 13 and 15 now)
fragment_lengths = [13, 15, 17]

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to create fragments and write to a FASTA file
def create_fasta(sequence, fraglength, protein_name, output_dir):
    start = 1  # Start position for the first fragment
    increment = 1  # Step increment for creating new fragments
    length = len(sequence)
    fragments = []

    # Create a list of fragments of the specified length, spaced by the increment
    for i in range(0, length - fraglength + 1, increment):
        fragments.append(sequence[i:i + fraglength])

    # Define the output file name
    output_file = os.path.join(output_dir, f"{protein_name}_fragments_{fraglength}aa.fasta")

    # Create a FASTA file with fragments interspersed with 10 U's
    with open(output_file, 'w') as f:
        for i, frag in enumerate(fragments):
            # Write the FASTA header
            f.write(f'>{protein_name} {start + increment * i}-{start + increment * i + fraglength - 1}\n')

            n = 4  # Number of repetitions of the fragment with UUUUUUUUUU groups
            for _ in range(n):
                f.write(frag + 'UUUUUUUUUU')
            # Write the fragment followed by a newline
            f.write(frag + '\n')

    # Print the output path for debugging
    print(f"FASTA file saved to: {os.path.abspath(output_file)}")

# Loop through each fragment length and create a FASTA file
for fraglength in fragment_lengths:
    create_fasta(sequence, fraglength, protein_name, output_dir)

