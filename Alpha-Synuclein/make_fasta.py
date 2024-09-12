
#%%

#edit these for sequence, start value, and desired increment
#sequence = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
import argparse

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Process amino acid sequence and fragment length.")

# Add arguments for sequence and fragment length
parser.add_argument('--sequence', type=str, required=True, help='Amino acid sequence to be processed.')
parser.add_argument('--fraglength', type=int, required=True, help='Fragment length.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
sequence = args.sequence
fraglength = args.fraglength

start = 1
increment = 1

end = start + len(sequence) - 1
length = len(sequence)
fragments = []

#create list of 15 amino acid fragments spaced by increment
for i in range(0 , length-fraglength+1, increment):
    #print(resid[i:i+15])
    fragments.append(sequence[i:i+fraglength])
print(fragments)

# Define output file name and path
output_file = 'alpha_synuclein_fragments.fasta'

#create fasta file with fragments + 10*U + fragments + 10*U ... 
with open(output_file, 'w') as f:
    for i, frag in enumerate(fragments):
        f.write(f'>alpha_synuclein {start+increment*i}-{start+increment*i+fraglength-1}\n')

        n = 4 #number of UUUUUUUUUU groups
        for _ in range(n):  # Replace n with the desired number of iterations
            f.write(frag + 'UUUUUUUUUU')
        f.write(frag + '\n')

# List output file for bugfixing
output_path = os.path.abspath(output_file)
print(f"FASTA file saved to: {output_path}")
# %%





#%%
import argparse
import os

# Create the ArgumentParser object
parser = argparse.ArgumentParser(description="Process amino acid sequence and fragment length.")

# Add arguments for sequence, fragment length, and protein name
parser.add_argument('--sequence', type=str, required=True, help='Amino acid sequence to be processed.')
parser.add_argument('--protein_name', type=str, required=True, help='Name of the protein.')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
sequence = args.sequence
protein_name = args.protein_name

# Fragment lengths to be processed
fragment_lengths = [13, 15, 17, 19]

# Define the output directory based on protein_name
output_dir = os.path.join(protein_name, "FASTA")
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
            for _ in range(n):  # Replace n with the desired number of iterations
                f.write(frag + 'UUUUUUUUUU')
            # Write the fragment followed by a newline
            f.write(frag + '\n')

    # Print the output path for debugging
    print(f"FASTA file saved to: {os.path.abspath(output_file)}")

# Loop through each fragment length and create a FASTA file
for fraglength in fragment_lengths:
    create_fasta(sequence, fraglength, protein_name, output_dir)
