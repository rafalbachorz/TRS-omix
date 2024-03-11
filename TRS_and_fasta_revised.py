#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from modules.pytrsomix import TRScalculator, SeqAnalyzer, SeqProcessor
import matplotlib.pyplot as plt
from multiprocessing import Pool
import argparse


# In[2]:


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

# Check if results file exists in the current working directory
results_directory = os.path.join(os.getcwd(), f"{directory_name}_results")

results_file_path = os.path.join(os.getcwd(), results_directory)


# Define the directory path for the results
results_directory = os.path.join(os.getcwd(), f"{directory_name}_results")
results_file = f"{directory_name}_results.csv"
# Create the directory
SeqProcessor.ensure_directory_exists(results_directory)

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
SeqProcessor.ensure_directory_exists(trs_output_dir)
csv_file_path = os.path.join(trs_output_dir, csv_file_name)

# Save combined_results to CSV file
combined_results.to_csv(csv_file_path, index=False)

print(f"Results saved to {csv_file_path}")


# In[3]:


combined_results = pd.read_csv(csv_file_path)


# In[4]:


# Use SeqProcessing methods to calculate maximum lengths for SEQ_L and SEQ_R based on tmin
l_chars_max, r_chars_max = SeqProcessor.calculate_sequence_lengths(tmin)

# Adjust the user input to ensure it is within the allowed range
l_chars = SeqProcessor.adjust_input_to_range(f"Enter the length of sequence to extract from the start (up to {l_chars_max}): ", l_chars_max)
r_chars = SeqProcessor.adjust_input_to_range(f"Enter the length of sequence to extract from the start (up to {r_chars_max}): ", r_chars_max)

# Use SeqProcessing to extract sequences based on the adjusted lengths
combined_results = SeqProcessor.extract_sequences(combined_results, l_chars, r_chars)

# Generate the new directory name based on user inputs
new_results_directory = f"{results_directory}_L{l_chars}_R{r_chars}"
new_results_directory_path = os.path.join(os.path.dirname(results_directory), new_results_directory)

print(f"Results directory is set to: {new_results_directory}")

# From this point, you can continue with further processing or analysis steps
# such as saving the modified DataFrame to a new file in the new results directory
# or performing additional operations as needed.
if not os.path.exists(new_results_directory_path):
    # Rename the existing results directory
    os.rename(results_directory, new_results_directory_path)
    print(f"The results directory has been renamed to: {new_results_directory_path}")
else:
    print(f"Directory {new_results_directory_path} already exists. Consider using a different name or removing the existing directory.")
results_directory = new_results_directory_path 


# In[ ]:





# In[5]:


SeqProcessor.set_user_email()

ncbi_ids = combined_results["GENOME"].unique().tolist()

organism_map = SeqProcessor.fetch_organism_names(ncbi_ids)

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


# In[6]:


output_path = os.path.join(trs_output_dir, 'combined_sequences_unique.fasta')

fasta_combined_file_path = output_path

SeqProcessor.rename_sequences(fasta_file_path, output_path)


# In[8]:


fasta_file_path = os.path.join(trs_output_dir, 'combined_sequences.fasta')

output_path = os.path.join(results_directory, "cd-hit_results")

output_file = os.path.join(output_path, "combined_sequences_unique_cdhit")
cdhit_path = SeqProcessor.find_file_by_name(file_name='cd-hit-est', auto=True)
cdhit_path = Path(cdhit_path)
cdhit_path = cdhit_path.parent
if not os.path.exists(output_path):
    os.makedirs(output_path)
results_directory = SeqProcessor.run_cdhit_new(cdhit_path, input_file=fasta_combined_file_path, output_file=output_file,
                                           results_directory=results_directory, sc = 1)

trs_output_dir = os.path.join(results_directory,"TRS_output")


# In[9]:


input_file = os.path.join(results_directory,"cd-hit_results","combined_sequences_unique_cdhit.clstr")
output_file = os.path.join(results_directory,"cd-hit_results","combined_sequences_clusters.txt")
SeqProcessor.extract_sequence_ids(input_file,output_file)


# In[10]:


cluster_to_clean = os.path.join(results_directory,"cd-hit_results","combined_sequences_clusters.txt")


# In[11]:


SeqProcessor.clean_sequence_ids(cluster_to_clean)
cluster_cleaned = cluster_to_clean


# In[12]:


fasta_ids_to_remove = SeqProcessor.read_fasta_ids(cluster_cleaned)

fasta_combined_file_path = os.path.join(trs_output_dir, 'combined_sequences_unique.fasta')

filtered_sequences_directory = os.path.join(results_directory,"filtered_sequences")

SeqProcessor.ensure_directory_exists(filtered_sequences_directory)

filtered_fasta = os.path.join(filtered_sequences_directory,"filtered_sequences_combined_unique.fasta")

SeqProcessor.filter_fasta_file(fasta_combined_file_path,filtered_fasta, fasta_ids_to_remove)


# In[13]:



fasta_ids_to_include = SeqProcessor.read_fasta_ids(cluster_cleaned)

sequences_in_clusters = os.path.join(filtered_sequences_directory, "cluster_sequences_combined_unique.fasta")

unique_path = os.path.join(trs_output_dir, 'combined_sequences_unique.fasta')

SeqProcessor.filter_fasta_file_clusters(unique_path, sequences_in_clusters, fasta_ids_to_include)

blast_out_directory = os.path.join(results_directory,"blast_output")

SeqProcessor.ensure_directory_exists(blast_out_directory)

print(f"""All results files are stored in {filtered_sequences_directory},
please continue analysis by blasting them with 100% identity in tabular format. Directory for the blast output
was created at {blast_out_directory}""")


# In[24]:


summary_dir = os.path.join(results_directory,"summary")
SeqProcessor.ensure_directory_exists(summary_dir)
summary_path = os.path.join(summary_dir,"summary.txt")
SeqProcessor.write_summary(summary_path,results_directory)


# In[16]:





# In[22]:





# In[23]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




