#!/usr/bin/env python
# coding: utf-8

# In[15]:


import cffi
import os
import numpy as np
import pandas as pd
from io import StringIO
from enum import Enum
from Bio import SeqIO
from Bio import Entrez
from pathlib import Path
from modules.auxiliary import DATA_SEQ_DIR
from modules.pytrsomix import TRScalculator, TRSanalyzer
import re


# '''
# TODO:
# Reuse functions some of them work on very similar data
# Refactor code before exporting to account for suboptimal inputs/do error handling
# Export functions to the modules - will probably need code rewrite
# Probably split code below into operations done on fasta files and on blast output - done
# (maybe provide example script for HPC clusters running on slurm) - pending
# Make alternative version that could run on Windows(no idea how to do that)
# Make functions accept args instead of input - currently on hold due to some kind of problem with storing args in variables
# Please try to break this code as much as possible in many cases i did not introduce any try/catch functionality - report them to me 
# Checkpointing functions to allow restarting from a given step - this is pretty hard with my current skillset 
# '''

# In[9]:


def ensure_directory_exists(directory_path):
    """
    Ensure that the specified directory exists. If it does not exist, it is created.

    Parameters:
    - directory_path: The path of the directory to check and potentially create.
    """
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory {directory_path} created.")


# In[11]:


# Hopefully better version of the above
# Prompt user for directory path
print("Please provide the directory path where your .fasta files are located.")
print("You can provide either a full path or a relative path.")

# Loop until a valid directory path is provided
while True:
    directory_path = input("Directory path: ")
    directory_path = os.path.abspath(directory_path)
    if os.path.exists(directory_path):
        break
    else:
        print("Directory does not exist. Please provide a valid directory path.")

while True:
    try:
        tmin = int(input("Enter the value for minimum length of sequence: "))
        tmax = int(input("Enter the value for maximum length of sequence: "))
        mode = int(input("Enter the value for mode to be used: "))
        break
    except ValueError:
        print("Invalid input. Please enter a valid integer.")


# Get directory name
directory_name = os.path.basename(directory_path)


# Define the directory path for the results
results_directory = os.path.join(os.getcwd(), f"{directory_name}_results")
results_file = f"{directory_name}_results.csv"
# Create the directory
os.makedirs(results_directory, exist_ok=True)

# Define the file path for the results within the directory
results_file_path = os.path.join(results_directory, results_file)

if os.path.exists(results_file_path):
    # Ask user if they want to continue
    choice = input(f"Results file '{results_file}' already exists in the current directory. Do you want to continue? (y/n): ")
    if choice.lower() != 'y':
        print("Exiting program.")
        exit()
        
# Get list of .fasta files in the directory
fasta_files = [file for file in os.listdir(directory_path) if file.endswith(".fasta")]

# Initialize a list to store TRScalculator instances
trs_calculators = []

# Iterate over each .fasta file
for fasta_file in fasta_files:
    # Get absolute path of the fasta file
    fasta_path = os.path.join(directory_path, fasta_file)

    # Check if file exists
    if not os.path.exists(fasta_path):
        print(f"File '{fasta_file}' does not exist. Skipping...")
        continue
    
    # Define the file path for trs.txt (assuming it's located in the working directory)
    trs_file = os.path.abspath("trs.txt").encode()
    
    # Create TRScalculator instance and calculate TRS
    trs_calculator = TRScalculator(sequence=fasta_path.encode(), trs=trs_file, tmin=tmin, tmax=tmax, mode=mode)
    trs_calculator.calculate()
    
    # Append the TRScalculator instance to the list
    trs_calculators.append(trs_calculator)

# Initialize an empty list to store the results
results_list = []

# Iterate over each TRScalculator instance
for trs_calculator in trs_calculators:
    # Extract the result from the calculator
    result = trs_calculator.Result
    
    # Append the result to the list
    results_list.append(result)

# Concatenate all results into a single DataFrame
combined_results = pd.concat(results_list, ignore_index=True)

#Remove ">" from >SEQ column
combined_results[">SEQ"] = combined_results[">SEQ"].str.replace(">", "")

# Display the combined results
combined_results

# Define the CSV file name and directory 
csv_file_name = directory_name + "_results.csv"
trs_output_dir = os.path.join(results_directory, "TRS_output")
os.makedirs(trs_output_dir, exist_ok=True)
csv_file_path = os.path.join(trs_output_dir, csv_file_name)

# Save combined_results to CSV file
combined_results.to_csv(csv_file_path, index=False)

print(f"Results saved to {csv_file_path}")


# In[16]:





# In[ ]:





# In[19]:


# Now we can load
combined_results = pd.read_csv(csv_file_path)
combined_results


# In[20]:


#Make get maximum values for user input based on what they gave as tmin
def calculate_sequence_lengths(tmin):
    l_chars_max = tmin // 2
    r_chars_max = tmin // 2
    return l_chars_max, r_chars_max

l_chars_max, r_chars_max = calculate_sequence_lengths(tmin)

print(f"The maximum available length for SEQ_L and SEQ_R is {l_chars_max} and {r_chars_max} respectively.")

l_chars = int(input(f"Enter the length of sequence to extract from the start (up to {l_chars_max}): "))
r_chars = int(input(f"Enter the length of sequence to extract from the end (up to {r_chars_max}): "))

# Not very elegant but catch all solution to the problem of bad inputs (i know it does not account for string will work on it)
def adjust_input_to_range(user_input, max_val):
    if user_input > max_val:
        print(f"Your input was adjusted to {max_val}.")
        return max_val
    return user_input

l_chars = adjust_input_to_range(l_chars, l_chars_max)
r_chars = adjust_input_to_range(r_chars, r_chars_max)

def extract_sequences(combined_results, l_chars, r_chars):
    # Create new columns SEQ_L and SEQ_R by slicing the sequence column
    combined_results['SEQ_L'] = combined_results['>SEQ'].str.slice(0, l_chars)
    combined_results['SEQ_R'] = combined_results['>SEQ'].str.slice(-r_chars)    
    return combined_results

combined_results = extract_sequences(combined_results, l_chars, r_chars)

# Generate the new directory name based on user inputs
new_results_directory = f"{results_directory}_L{l_chars}_R{r_chars}"

# Full path for the new directory
new_results_directory_path = os.path.join(os.path.dirname(results_directory), new_results_directory)

# Check if the new directory name already exists to avoid overwriting
if not os.path.exists(new_results_directory_path):
    # Rename the existing results directory
    os.rename(results_directory, new_results_directory_path)
    print(f"The results directory has been renamed to: {new_results_directory_path}")
else:
    print(f"Directory {new_results_directory_path} already exists. Consider using a different name or removing the existing directory.")
results_directory = new_results_directory_path    


# In[21]:


# Function to validate email format
def validate_email(email):
    """
    Validate email format.

    Args:
        email (str): Email address to validate.

    Returns:
        bool: True if the email format is valid, False otherwise.
    """
    if re.match(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$", email):
        return True
    else:
        return False

# Prompt the user to input their email address
while True:
    print("To further protect your email address, consider using a temporary or disposable email address.")
    user_email = input("Enter your email address to access Entrez: ")
    # Validate email format
    if validate_email(user_email):
        break
    else:
        print("Invalid email address format. Please enter a valid email address.")

# Set the user's email address (required by NCBI)
Entrez.email = user_email

# Function to fetch organism names from NCBI based on NCBI IDs
def fetch_organism_names(ncbi_ids):
    """
    Fetch organism names from NCBI based on NCBI IDs.

    Args:
        ncbi_ids (list): List of NCBI IDs.

    Returns:
        dict: Dictionary mapping NCBI IDs to organism names.
    """
    organism_map = {}  # Initialize a dictionary to store the organism names
    for ncbi_id in ncbi_ids:
        try:
            # Fetch the record from NCBI
            handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
            record = handle.read()
            # Extract the organism information from the record
            for line in record.splitlines():
                if line.startswith("  ORGANISM"):
                    organism = line.split("ORGANISM")[1].strip()
                    # Replace spaces with underscores
                    organism = organism.replace(" ", "_")
                    organism_map[ncbi_id] = organism
                    break  # Once organism information is found, break the loop
        except Exception as e:
            print(f"Error fetching organism name for NCBI ID {ncbi_id}: {e}")
    return organism_map

# Read the DataFrame with the combined_results
ncbi_ids = combined_results["GENOME"].unique().tolist()

# Fetch organism names from NCBI
organism_map = fetch_organism_names(ncbi_ids)


# In[22]:


# Map NCBI IDs to taxonomic names
combined_results['Taxonomic_Name'] = combined_results['GENOME'].map(organism_map)

# Identify unmatched genomes
unmatched_genomes = combined_results[combined_results['Taxonomic_Name'].isnull()]['GENOME'].unique()

if len(unmatched_genomes) > 0:
    print(f"Warning: Some genome IDs could not be matched with taxonomic names: {unmatched_genomes}")

# Extract sequences and create sequence IDs
combined_results['L_id'] = combined_results['Taxonomic_Name'] + '_L' + combined_results['L-No'].astype(str)
combined_results['R_id'] = combined_results['Taxonomic_Name'] + '_R' + combined_results['R-No'].astype(str)

# Select relevant columns
sequences_df = combined_results[['SEQ_L', 'SEQ_R', 'L_id', 'R_id']]

# Save both left and right sequences to a single FASTA file will help with cd-hit
trs_output_dir = os.path.join(results_directory,"TRS_output")
fasta_file_path = os.path.join(trs_output_dir, 'combined_sequences.fasta')
with open(fasta_file_path, 'w') as fasta_file:
    for _, row in sequences_df.iterrows():
        # Write left sequence
        fasta_file.write(f'>{row["L_id"]}\n')
        fasta_file.write(f'{row["SEQ_L"]}\n')
        # Write right sequence
        fasta_file.write(f'>{row["R_id"]}\n')
        fasta_file.write(f'{row["SEQ_R"]}\n')


# In[23]:


#There is possibility of several sequences having the same name as such we need to make sure that all ids are unique
def rename_sequences(input_file, output_file):
    """
    Renames sequences in a FASTA file such that each pair of sequences (L and R) shares the same index.

    Parameters:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file with renamed sequences.
    """
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        # Initialize a counter for pairs
        pair_index = 1
        # Process sequences in pairs
        for record in SeqIO.parse(input_handle, "fasta"):
            original_id = record.id
            # Determine whether the sequence is part of a pair (L/R) and assign the same index
            # Assuming sequences come in consecutive pairs (L followed by R)
            new_id = f"{original_id}_{pair_index}"
            record.id = new_id
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")
            
            # Increment the pair index after every second sequence to keep the index the same for pairs
            if 'R' in original_id:
                pair_index += 1

# Call the function to rename sequences
output_path = os.path.join(trs_output_dir, 'combined_sequences_unique.fasta')   
rename_sequences(fasta_file_path, output_path)
fasta_combined_file_path = output_path
# This function incorporates index into sequence name making it easy to find the original sequence


# In[25]:


def run_cdhit(cdhit_path, input_file, output_file, c=None, d=0, m="t", g=0, G=1, sc=0, results_directory="results"):
    """
    Enhanced run CD-HIT program from Python, automatically adjusting 'n' based on sequence identity threshold 'c'.

    Parameters:
        cdhit_path (str): Path to the CD-HIT executable. If None, will use a globally stored path or prompt the user.
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output file.
        c (float): Sequence identity threshold (between 0.4 and 1.0). If None, the user will be prompted.
        d (int): Bandwidth of alignment, default is 0.
        m (str): Memory limit, default is "t" for unlimited.
        g (int), G (int), sc (int): Additional CD-HIT parameters with their default values.
        results_directory (str): Directory to store the results. Default is "results".
    """
    import os
    import subprocess

    global cdhit_path_global

    # Function to adjust 'n' based on 'c'
    def adjust_word_length(c):
        thresholds_n_values = [
            (0.95, 10),
            (0.90, 8),
            (0.88, 7),
            (0.85, 6),
            (0.80, 5),
            (0.75, 4),
        ]
        if c == 1.0:
            return 11
        for lower_bound, n_value in thresholds_n_values:
            if c >= lower_bound:
                return n_value
        return None  # Fallback, should not be reached if 'c' is within specified bounds

    # Prompt for sequence identity threshold if not provided
    if c is None:
        c = float(input("Enter the sequence identity threshold (between 0.75 and 1.0): "))
        while not 0.75 <= c <= 1.0:
            print("Error: The sequence identity threshold must be between 0.75 and 1.0.")
            c = float(input("Enter the sequence identity threshold (between 0.75 and 1.0): "))

    # Dynamically adjust 'n' based on 'c'
    n = adjust_word_length(c)

    # Validate and use the CD-HIT path
    if cdhit_path is None and 'cdhit_path_global' in globals() and os.path.exists(os.path.join(cdhit_path_global, "cd-hit")):
        cdhit_path = cdhit_path_global
    elif cdhit_path is None or not os.path.exists(os.path.join(cdhit_path, "cd-hit")):
        cdhit_path = input("Enter the directory where CD-HIT is located: ")
        if not os.path.exists(os.path.join(cdhit_path, "cd-hit")):
            raise FileNotFoundError("CD-HIT executable not found in the specified directory.")
    cdhit_path_global = cdhit_path

    # Ensure the results directory exists
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    # Construct the CD-HIT command
    cmd = [
        os.path.join(cdhit_path, "cd-hit-est"),
        "-i", input_file,
        "-o", output_file,
        "-c", str(c),
        "-n", str(n),
        "-d", str(d),
        "-M", m,
        "-g", str(g),
        "-G", str(G),
        "-sc", str(sc)
    ]

    # Execute CD-HIT command and handle errors
    try:
        subprocess.run(cmd, check=True)
        new_results_directory = f"{results_directory}_c{c}"
        if not os.path.exists(new_results_directory):
            os.rename(results_directory, new_results_directory)
            results_directory = new_results_directory
            print(f"The results directory has been renamed to: {new_results_directory}")
        else:
            print(f"Warning: The directory {new_results_directory} already exists. Results directory was not renamed.")
    except subprocess.CalledProcessError as e:
        print(f"CD-HIT command failed: {e}")
        raise

    return results_directory

output_path = os.path.join(results_directory, "cd-hit_results")
output_file = os.path.join(output_path, "combined_sequences_unique_cdhit")

# Ensure the output directory exists
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Run CD-HIT and update results_directory with the potentially new directory name
results_directory = run_cdhit(cdhit_path=None, input_file=fasta_combined_file_path, output_file=output_file, results_directory=results_directory)
# results_directory now points to the updated directory name, if it was changed


# In[26]:


trs_output_dir = os.path.join(results_directory,"TRS_output")


# In[27]:


def extract_sequence_ids(input_file, output_file):
    """
    Extract sequence IDs from clusters with more than 2 sequences and save them to a file.

    Parameters:
        input_file (str): Path to the input cluster data file.
        output_file (str): Path to the output file to save sequence IDs.

    Returns:
        None
    """
    with open(input_file, "r") as f, open(output_file, "w") as outfile:
        current_cluster = []
        for line in f:
            if line.startswith(">Cluster"):
                if len(current_cluster) >= 2:
                    for seq in current_cluster:
                        outfile.write(seq)
                current_cluster = []  # Clear current cluster
            elif line.strip():  # Proceed if the line is not empty
                current_cluster.append(line)  # Add sequence to current cluster
        if len(current_cluster) >= 2:
            for seq in current_cluster:
                outfile.write(seq)
# Example usage
input_file = os.path.join(results_directory,"cd-hit_results","combined_sequences_unique_cdhit.clstr")
output_file = os.path.join(results_directory,"cd-hit_results","combined_sequences_clusters.txt")
extract_sequence_ids(input_file, output_file)


# In[28]:


def clean_sequence_ids(input_file):
    """
    Clean sequence IDs file by removing prefixes and suffixes.

    Parameters:
        input_file (str): Path to the input file with sequence IDs.

    Returns:
        None
    """
    with open(input_file, "r") as infile:
        lines = infile.readlines()

    with open(input_file, "w") as outfile:
        for line in lines:
            sequence_id = line.split(">")[1].split("...")[0].strip()
            outfile.write(">" + sequence_id + "\n")

# Example usage
input_file = output_file
clean_sequence_ids(input_file)


# In[29]:


"""
Remove entries in L/R_clusters from the .fasta file 
"""
def read_fasta_ids(filename):
    """
    Read and return a set of FASTA IDs from a given file, processing the file line by line
    to efficiently handle large files.

    Parameters:
    - filename: Path to the FASTA file.

    Returns:
    - A set of FASTA IDs found in the file.
    """
    fasta_ids = set()
    
    def process_line(line): #might be useful later
        """
        Process a single line from the FASTA file, adding the ID to the set if the line is a header.
        """
        if line.startswith('>'):
            fasta_id = line.strip()[1:]
            fasta_ids.add(fasta_id)

    try:
        with open(filename, 'r') as f:
            for line in f:
                process_line(line)
    except FileNotFoundError:
        print(f"Error: File {filename} not found. Do not move files during execution!")

    return fasta_ids


def filter_fasta_file(input_file, output_file, fasta_ids_to_remove, chunk_size=1024*1024):
    """
    Filter entries from a FASTA file in chunks, removing sequences with IDs in the given set.

    Parameters:
    - input_file: Path to the input FASTA file.
    - output_file: Path where the filtered FASTA file will be saved.
    - fasta_ids_to_remove: A set of FASTA IDs to be removed from the input file.
    - chunk_size: Size of the chunk to read at a time (in bytes). Default is 1MB.
    """
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            fasta_entry = []  # Stores lines of a single FASTA entry
            write_entry = True  # Determines whether the current entry should be written to the output

            for line in f_in:
                if line.startswith('>'):  # Start of a new FASTA entry
                    if fasta_entry:  # If there's an entry to process
                        if write_entry:  # If the previous entry is not in the removal list, write it to the file
                            f_out.writelines(fasta_entry)
                        fasta_entry = []  # Reset for the next entry
                        write_entry = True  # Reset flag

                    current_fasta_id = line.strip()[1:]
                    if current_fasta_id in fasta_ids_to_remove:
                        write_entry = False  # Mark for non-writing if ID is in the removal list

                fasta_entry.append(line)  # Add current line to the entry

            if fasta_entry and write_entry:  # Process the last entry if there's one
                f_out.writelines(fasta_entry)

    except FileNotFoundError:
        print(f"Error: File {input_file} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Read fasta IDs from the L_clusters file
fasta_ids_to_remove = read_fasta_ids(input_file)
fasta_combined_file_path = os.path.join(trs_output_dir, 'combined_sequences_unique.fasta') 

#Make sure the output directory exists
filtered_sequences_directory = os.path.join(results_directory,"filtered_sequences")
ensure_directory_exists(filtered_sequences_directory)

# Filter the sequences_L_unique.fasta file
filtered_fasta = os.path.join(filtered_sequences_directory,"filtered_sequences_combined_unique.fasta")
filter_fasta_file(fasta_combined_file_path,filtered_fasta, fasta_ids_to_remove)

print(f"Filtered file created successfully at {filtered_fasta}!")


# In[30]:


def filter_fasta_file_clusters(input_file, output_file, fasta_ids_to_include):
    """
    Filters a FASTA file to only include sequences with IDs present in the given set.

    Parameters:
    - input_file: Path to the input FASTA file.
    - output_file: Path where the filtered FASTA file will be saved.
    - fasta_ids_to_include: A set of FASTA IDs that should be included in the output file.
    """
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            include_entry = False
            for line in f_in:
                if line.startswith('>'):
                    current_fasta_id = line.strip()[1:]
                    if current_fasta_id in fasta_ids_to_include:
                        include_entry = True
                        f_out.write(line)  # Include the header line for matching IDs
                    else:
                        include_entry = False  # Do not include entries for non-matching IDs
                elif include_entry:
                    f_out.write(line)  # Include the sequence lines for matching IDs
    except FileNotFoundError:
        print(f"Error: The file {input_file} was not found. Please check the file path.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# In[31]:


# Ensure the filtered_sequences directory exists
ensure_directory_exists(filtered_sequences_directory)

# Assuming the correct file paths are defined here for cluster files
clusters_file = os.path.join(results_directory,"cd-hit_results","combined_sequences_clusters.txt") 
# Read fasta IDs from the clusters file
fasta_ids_to_include = read_fasta_ids(clusters_file)

# Define the output paths and paths of initial sequences
sequences_in_clusters = os.path.join(filtered_sequences_directory, "cluster_sequences_combined_unique.fasta")
unique_path = os.path.join(trs_output_dir, 'combined_sequences_unique.fasta')


# Filter the combined_sequences_unique.fasta file
filter_fasta_file_clusters(unique_path, sequences_in_clusters, fasta_ids_to_include)
print(f"Cluster sequences file created successfully! {sequences_in_clusters}")

blast_out_directory = os.path.join(results_directory,"blast_output")
ensure_directory_exists(blast_out_directory)

print(f"""All results files are stored in {filtered_sequences_directory},
please continue analysis by blasting them with 100% identity in tabular format. Directory for the blast output
was created at {blast_out_directory}""")


# In[ ]:




