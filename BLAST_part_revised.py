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
from modules.pytrsomix import TRScalculator, SeqAnalyzer , SeqProcessor
import matplotlib.pyplot as plt
from multiprocessing import Pool
import glob
import re
import fnmatch


# In[2]:


#Look for blast_output directory
blast_output_path = SeqProcessor.find_directory_by_name('blast_output')
blast_output_path= Path(blast_output_path)
results_directory = blast_output_path.parent
if blast_output_path:
    print(f"Selected directory: {blast_output_path}")
else:
    # If blast_output_path is either None or an empty string, prompt the user for the path.
    print("No 'blast_output' directory selected or found.")
    blast_output_path = input("Please provide the path to the 'blast_output' directory: ")

#Renames files to .txt and find unique sequence ID - accession pairs  
SeqProcessor.process_blast_files_in_directory(blast_output_path)


# In[3]:


'''
Finds the taxonomy database 
on the users computer name is hardcoded initially as this is the name of the unpacked database from NCBI's FTP
'''
taxonomy_file = SeqProcessor.find_file_by_name('nucl_gb.accession2taxid',auto=True)
if taxonomy_file:
    print(f"Selected file: {taxonomy_file}.")
else:
    print(f"No taxonomy selected or found.")
    taxonomy_file = input("Please provide the path to nucl_gb.accession2taxid or equivalent: ")
'''
Get accessions present in blast files and filter the taxonomy_file using them this way we have all accession - taxid
Pairs in out dataset and there is not need to access such a large database.
In addition as you can see I'm using chunks to process the initial database, that size could be altered but i found 
50000 to be sweet spot. Loading progress bar was added for this operation as well.
'''
accessions = SeqProcessor.collect_accessions_from_blast_files(blast_output_path)
tax_df = SeqProcessor.filter_taxonomy_file(taxonomy_file,accessions,50000)


'''
Saving the filtered database to a csv file for use later - NOT IMPLEMENTED
'''
output_file = os.path.join(blast_output_path, "taxonomy_filtered.csv")
tax_df.to_csv(output_file, index=False)
print(f"tax_df saved to {output_file}")


# In[4]:


'''
Quite self-explanatory makes sure that each taxid - accession pair gets put into dictionary 
only once i make use of sets here to decrease it's size and repeated information
'''
taxid_accessions_dict = {}
for index, row in tax_df.iterrows():
    accession_column = tax_df.columns[0]  # Extract accession column dynamically
    taxid_column = tax_df.columns[1]  # Extract taxid column dynamically
    
    accession = row[accession_column]
    taxid = row[taxid_column]
    
    if taxid in taxid_accessions_dict:
        taxid_accessions_dict[taxid].append(accession)
    else:
        taxid_accessions_dict[taxid] = [accession]


# In[5]:


# Ensure the directory for modified files exists
modified_blast_path = os.path.join(blast_output_path, "modified_blast")
SeqProcessor.ensure_directory_exists(modified_blast_path)  # Ensure this function creates the directory if it doesn't exist

# Optimized: Invert the taxid_accessions_dict for efficient lookup
'''
This step assumes each accession maps to exactly one taxid which is true (but that matching is not always accurate)
this drawback is attributed to nucleotide database structure and is not something that should alter the results in any
significant way as long as we remember to add the top-level species to the dictionary later
'''
accession_to_taxid = {accession: taxid for taxid, accessions in taxid_accessions_dict.items() for accession in accessions}

def map_accession_to_taxid(accession):
    """Map accession to TaxID using the optimized lookup dictionary."""
    return accession_to_taxid.get(accession, '')  # Return empty string if accession not found

# Iterate over files in the input directory
for filename in os.listdir(blast_output_path):
    if filename.endswith(".txt"):
        input_file_path = os.path.join(blast_output_path, filename)
        
        # Read the file into a DataFrame
        df = pd.read_csv(input_file_path, sep='\t', header=None)  # No column names specified
        
        # Map accessions to TaxIDs and add them as a new column after the last existing column
        df[df.shape[1]] = df[1].map(map_accession_to_taxid)  # Assuming column 1 contains the accessions
        
        # Construct output file path
        modified_file_path = os.path.join(modified_blast_path, f"taxids_{filename}")
        
        # Save the modified DataFrame to a new file in the output directory
        df.to_csv(modified_file_path, sep='\t', index=False, header=False)


# In[6]:


'''This is still experimental but should allow to load interiors.txt equivalent into a dataframe to construct the
dictionary of genomes which we analyze.'''
print("Searching for the results file....")
results_file = SeqProcessor.find_file_by_name('*_results.csv',folder = results_directory, auto=True)
if results_file:
    print(f"Selected file: {results_file}")
else:
    # If results_file is either None or an empty string, prompt the user for the path.
    print("No '*_results.csv' file selected or found.")
    blast_output_path = input("Please provide the path to the current experiments result file: ")

combined_results = pd.read_csv(results_file)


# In[7]:


'''
Set mail for Entrez, extract unique genome ids from dataframe
'''
SeqProcessor.set_user_email()
ncbi_ids = combined_results["GENOME"].unique().tolist()
'''
Retrieves organism names first, then retrieves taxids for those IDs, creates a new dict name:taxid
'''
tax_map = SeqProcessor.fetch_organism_taxids(ncbi_ids)


# In[8]:


'''
Filtering functions gets rid of _parts_ that contain numbers or the capital letter(save for the first one)
Should get rid of strains and make it possible to find the taxid of the species in the next step
'''
filtered_organism_taxid_map = {SeqProcessor.filter_key_parts(key): value for key, value in tax_map.items()}


# In[9]:


filtered_organism_taxid_map 


# In[10]:


species_info = SeqProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)


# In[11]:


SeqProcessor.interact_and_update_dict(filtered_organism_taxid_map)


# In[12]:


species_info = SeqProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)
species_info


# In[13]:


species_info = {SeqProcessor.filter_key_parts(key): value for key, value in species_info.items()}


# In[14]:


modified_blast_path = os.path.join(blast_output_path, "modified_blast")
nan_file = os.path.join(modified_blast_path,"NaN acessions.csv")
# Construct dictionary of all taxid - acessions pairs in our data
results_dict = SeqProcessor.construct_dict_from_files(modified_blast_path,nan_file)


# In[ ]:





# In[15]:


exceptions = SeqProcessor.ask_for_exceptions()
print(f"Exceptions to filtering added: {exceptions}")


# In[16]:


for key, value_set in results_dict.items():
    results_dict[key] = {int(val) for val in value_set}    
    
for key, value_set in species_info.items():
    species_info[key] = {int(val) for val in value_set}    


# In[17]:


exceptions


# In[18]:


filtered_keys = SeqProcessor.filter_with_exceptions(results_dict,species_info,exceptions)

filtered_keys


# In[19]:


filtered_keys_final = SeqProcessor.unpack_single_element_sets(filtered_keys)
filtered_keys_final


# In[20]:


SeqProcessor.process_files_with_filter(modified_blast_path,filtered_keys_final)


# In[21]:


summary_dir = os.path.join(results_directory,"summary")
SeqProcessor.ensure_directory_exists(summary_dir)
summary_path = os.path.join(summary_dir,"summary.txt")
SeqProcessor.write_summary(summary_path,results_directory)


# In[22]:


fasta_processed_path = SeqProcessor.find_file_by_name('unique_taxids_filtered_sequences_combined_unique_blastn_out.txt',auto= True, folder= results_directory)
fasta_processed_path = Path(fasta_processed_path)
fasta_processed_path = fasta_processed_path.parent
fasta_processed_path


# In[24]:


SeqProcessor.filter_and_overwrite_files_in_directory(fasta_processed_path,file_pattern="*.txt")


# In[28]:


filtered_fasta_path = SeqProcessor.find_directory_by_name_new('filtered_sequences', auto= True, folder = results_directory)


# In[34]:


filtered_fasta_file = SeqProcessor.find_file_by_name('filtered_sequences_combined_unique.fasta', auto=True, folder = results_directory)


# In[25]:


unique_directory = SeqProcessor.find_directory_by_name_new('unique_sequences', auto = True, folder=results_directory)


# In[27]:


sequence_ids_dict, special_dict = SeqProcessor.read_sequence_ids(unique_directory)


# In[37]:





# In[38]:


output_directory = os.path.join(results_directory,"final_output")
SeqProcessor.ensure_directory_exists(output_directory)
SeqProcessor.filter_fasta_file(filtered_fasta_file,sequence_ids_dict,special_dict,output_directory)


# In[ ]:




