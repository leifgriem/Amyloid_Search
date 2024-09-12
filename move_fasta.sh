#!/bin/bash

# Define the search pattern and destination directory
SEARCH_DIR="./"
REMOTE_DIR="pod:/home/leifgriem/Amyloid_Project/FASTA"

# Debugging: List files being searched
echo "Searching for files in $SEARCH_DIR with patterns *13aa*.fasta, *15aa*.fasta, and *17aa*.fasta"
find "$SEARCH_DIR" -type f \( -name "*13aa*.fasta" -o -name "*15aa*.fasta" -o -name "*17aa*.fasta" \) -print

# Find and move .fasta files containing "13aa" or "15aa" or "17aa" in their name
find "$SEARCH_DIR" -type f \( -name "*13aa*.fasta" -o -name "*15aa*.fasta" -o -name "*17aa*.fasta" \) | while read file; do
    # Extract the filename from the full path
    FILENAME=$(basename "$file")

    # Debug print to check if it's processing each file
    echo "Checking $FILENAME"

    # Check if the file already exists in the remote directory
    if ssh pod "[ -f /home/leifgriem/Amyloid_Project/FASTA/$FILENAME ]" 2>/dev/null; then
        echo "$FILENAME already exists in the remote directory, skipping transfer."
    else
        echo "Transferring $FILENAME to $REMOTE_DIR"
        scp "$file" "$REMOTE_DIR"
    fi
done

echo "File transfer completed."

