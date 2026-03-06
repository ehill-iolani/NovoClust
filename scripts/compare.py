from collections import defaultdict
from Bio import SeqIO   # requires biopython
from Bio import AlignIO   # requires biopython
import sys
import os

# Get data paths
# os.listdir(data_dir)

# Add option to specify a different directory for the input data
if len(sys.argv) > 1:
    data_dir = sys.argv[1]
    data = os.listdir(data_dir)
else:
    data_dir = "cluster_analysis/variants_over1pct"
    data = os.listdir(data_dir)
    
# Add option to specify a different directory for the output data
if len(sys.argv) > 2:
    output_dir = sys.argv[2]
else:
    output_dir = "compare_analysis"

# Only find "consolidated.fasta" files
data = [file for file in data if file.endswith("consolidated.fasta")]

file_names = []
for file in data:
    if file.endswith(".fasta") or file.endswith(".fa"):
        file_names.append(file)


# Function to rename FASTA headers by combining headers of identical sequences and adding the basename of the input file
def rename_fasta_headers(input_fasta, output_fasta):
    # Create a dictionary to store sequences and their corresponding headers
    seq_dict = defaultdict(list)

    # Read the input FASTA file and populate the dictionary
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_dict[str(record.seq)].append(record.id)

    # add basename of the input file to the headers
    for seq in seq_dict:
        input_fasta_basename = os.path.splitext(os.path.basename(input_fasta))[0]
        seq_dict[seq] = [f"{input_fasta_basename}_{header}" for header in seq_dict[seq]]

    # Write the output FASTA file with renamed headers
    with open(output_fasta, "w") as output_handle:
        for seq, headers in seq_dict.items():
            new_header = "_".join(headers)  # Combine headers of identical sequences
            output_handle.write(f">{new_header}\n{seq}\n")
            
# Make output directories if they don't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)        

if not os.path.exists(os.path.join(output_dir, "cleaned_data")):
    os.makedirs(os.path.join(output_dir, "cleaned_data"))

if not os.path.exists(os.path.join(output_dir, "aligned_data")):
    os.makedirs(os.path.join(output_dir, "aligned_data"))

if not os.path.exists(os.path.join(output_dir, "clustered_data")):
    os.makedirs(os.path.join(output_dir, "clustered_data"))
    
if not os.path.exists(os.path.join(output_dir, "clustered_data/seq_clusters")):
    os.makedirs(os.path.join(output_dir, "clustered_data/seq_clusters"))

if not os.path.exists(os.path.join(output_dir, "combined_cleaned_data")):
    os.makedirs(os.path.join(output_dir, "combined_cleaned_data"))

# Loop through the list of file names and rename headers for each file
for file_name in file_names:
    input_fasta = os.path.join(data_dir, file_name)
    output_fasta = os.path.join(output_dir, "cleaned_data", file_name)
    rename_fasta_headers(input_fasta, output_fasta)

# Combine all the cleaned fasta files into one file
with open(os.path.join(output_dir, "combined_cleaned_data", "combined.fasta"), "w") as combined_handle:
    for file_name in file_names:
        cleaned_fasta = os.path.join(output_dir, "cleaned_data", file_name)
        with open(cleaned_fasta, "r") as cleaned_handle:
            combined_handle.write(cleaned_handle.read())

# Perform MAFFT alignment on the combined fasta file    
os.system(f"mafft --preservecase --auto {os.path.join(output_dir, 'combined_cleaned_data', 'combined.fasta')} > {os.path.join(output_dir, 'aligned_data', 'aligned.fasta')}")

# end trim the alignment
os.system(f"python scripts/endtrim_alignment.py {os.path.join(output_dir, 'aligned_data', 'aligned.fasta')} {os.path.join(output_dir, 'aligned_data', 'aligned_trimmed.fasta')}")

# ungap the trimmed alignment
os.system(f"python scripts/ungap_fasta.py {os.path.join(output_dir, 'aligned_data', 'aligned_trimmed.fasta')} {os.path.join(output_dir, 'aligned_data', 'aligned_trimmed_ungapped.fasta')}")

# Perform clustering on the trimmed alignment
os.system(f"vsearch --uc {os.path.join(output_dir, 'clustered_data', 'derep.uc')} --derep_fulllength {os.path.join(output_dir, 'aligned_data', 'aligned_trimmed_ungapped.fasta')} --output {os.path.join(output_dir, 'clustered_data', 'derep.fasta')} --sizeout")

# Read in the UC file
clusters = defaultdict(list)
with open(os.path.join(output_dir, "clustered_data", "derep.uc"), "r") as uc_file:
    for line in uc_file:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        cluster_id = parts[1]
        sequence_id = parts[8]
        clusters[cluster_id].append(sequence_id)    
        
# Find clusters with > 1 sequence
clusters_with_multiple_sequences = {cluster_id: seq_ids for cluster_id, seq_ids in clusters.items() if len(seq_ids) > 1}

# Remove the duplicate sequence IDs in dict elements
for cluster_id, seq_ids in clusters_with_multiple_sequences.items():
    clusters_with_multiple_sequences[cluster_id] = list(set(seq_ids))

# Print the clusters with > 1 sequence to one file
with open(os.path.join(output_dir, "clustered_data", "clusters_with_multiple_sequences.txt"), "w") as output_file:
    for cluster_id, seq_ids in clusters_with_multiple_sequences.items():
        output_file.write(f"Cluster {cluster_id}:\n")
        for seq_id in seq_ids:
            output_file.write(f"  {seq_id}\n")
            
# Write a new file for each cluster of sequences with > 1 sequence, containing the sequences in that cluster
for cluster_id, seq_ids in clusters_with_multiple_sequences.items():
    with open(os.path.join(output_dir, "clustered_data/seq_clusters", f"cluster_{cluster_id}.fasta"), "w") as output_file:
        for seq_id in seq_ids:
            # Find the sequence in the aligned fasta file
            for record in SeqIO.parse(os.path.join(output_dir, "combined_cleaned_data", "combined.fasta"), "fasta"):
                if record.id == seq_id:
                    SeqIO.write(record, output_file, "fasta")
                    break
