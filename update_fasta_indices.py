# Update FASTA Indices
import re
import argparse
import tempfile
import os

def update_fasta_indices(input_fasta, increment):
    # Create a temporary file
    with tempfile.NamedTemporaryFile('w', delete=False) as temp_file:
        temp_file_name = temp_file.name
        
        # Open the input FASTA file for reading
        with open(input_fasta, 'r') as infile:
            for line in infile:
                if line.startswith('>'):
                    # Match indices in the header (e.g., '>SAA 1-15')
                    match = re.search(r'(\d+)-(\d+)', line)
                    if match:
                        # Increment both indices by the specified increment
                        start = int(match.group(1)) + increment
                        end = int(match.group(2)) + increment
                        # Replace the old indices with the new ones
                        new_header = re.sub(r'\d+-\d+', f'{start}-{end}', line)
                        temp_file.write(new_header)
                    else:
                        temp_file.write(line)
                else:
                    # Write sequence lines as they are
                    temp_file.write(line)

    # Overwrite the original input file with the contents of the temp file
    os.replace(temp_file_name, input_fasta)

# Command-line interface using argparse
parser = argparse.ArgumentParser(description='Update the indices in a FASTA file header by a specified increment.')
parser.add_argument('--input', required=True, help='Path to the input FASTA file')
parser.add_argument('--increment', type=int, default=1, help='The increment by which to adjust the indices (default is 1)')

args = parser.parse_args()

# Run the function to update indices with the specified increment
update_fasta_indices(args.input, args.increment)

print(f"FASTA file indices updated by {args.increment} and the input file has been overwritten.")

