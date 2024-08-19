
#%%

#edit these for sequence, start value, and desired increment
sequence = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
start = 1
increment = 1
fraglength = 15

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

