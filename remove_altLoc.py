import os
import tempfile
import shutil  
import argparse  

def process_pdb_file_in_place(file_path):
    """
    Process a PDB file to remove all altLocs except 'A' and save the changes in the same file.
    """
    # Create a temporary file to write the changes
    with tempfile.NamedTemporaryFile('w', delete=False) as temp_file:
        temp_file_path = temp_file.name

        with open(file_path, 'r') as infile:
            for line in infile:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # AltLoc is located at position 17 in PDB ATOM lines (1-based index)
                    altLoc = line[16]  # 0-based index, 16 corresponds to position 17
                    # Only keep lines where altLoc is 'A' or blank (primary conformation)
                    if altLoc == 'A' or altLoc == ' ':
                        temp_file.write(line)
                else:
                    # Write non-ATOM/HETATM lines as-is
                    temp_file.write(line)

    # Use shutil.move instead of os.replace to ensure cross-device compatibility
    shutil.move(temp_file_path, file_path)
    print(f"Processed: {file_path}")

def process_directory_in_place(input_dir):
    """
    Process all PDB files in the input directory in place (overwrite the same files).
    """
    for filename in os.listdir(input_dir):
        if filename.endswith(".pdb"):
            file_path = os.path.join(input_dir, filename)
            process_pdb_file_in_place(file_path)

# Set up argparse for command-line argument parsing
parser = argparse.ArgumentParser(description="Remove altLocs except 'A' from PDB files in a directory.")
parser.add_argument('--input_dir', type=str, required=True, help="Path to the directory containing PDB files.")

# Parse the arguments
args = parser.parse_args()

# Run the altLoc removal process in place
process_directory_in_place(args.input_dir)
