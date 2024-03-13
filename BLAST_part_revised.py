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


# In[4]:


#Look for blast_output directory
blast_output_path = SeqProcessor.find_directory_by_name('blast_output')


if blast_output_path:
    print(f"Selected directory: {blast_output_path}")
else:
    # If blast_output_path is either None or an empty string, prompt the user for the path.
    print("No 'blast_output' directory selected or found.")
    blast_output_path = input("Please provide the path to the 'blast_output' directory: ")

#Renames files to .txt and find unique sequence ID - accession pairs  
SeqProcessor.process_blast_files_in_directory(blast_output_path)

blast_output_path= Path(blast_output_path)
results_directory = blast_output_path.parent


# In[ ]:


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


# In[ ]:


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


# In[ ]:


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


# In[ ]:


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


# In[ ]:


'''
Set mail for Entrez, extract unique genome ids from dataframe
'''
SeqProcessor.set_user_email()
ncbi_ids = combined_results["GENOME"].unique().tolist()
'''
Retrieves organism names first, then retrieves taxids for those IDs, creates a new dict name:taxid
'''
tax_map = SeqProcessor.fetch_organism_taxids(ncbi_ids)


# In[ ]:


'''
Filtering functions gets rid of _parts_ that contain numbers or the capital letter(save for the first one)
Should get rid of strains and make it possible to find the taxid of the species in the next step
'''
filtered_organism_taxid_map = {SeqProcessor.filter_key_parts(key): value for key, value in tax_map.items()}


# In[ ]:


filtered_organism_taxid_map 


# In[ ]:


'''Adds additional taxids after clearing out the strain associated info in current klebsiella dataset 
because of that Klebsiella_pneumoniae_subsp._pneumoniae has two taxids one strain specific and the one associated with subspecies'''
species_info = SeqProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)


# In[ ]:


'''Additional NCBI IDs can be provided to include their taxids in the dictionary used for filtering down the line
Most useful to add top level species in case of klebsiella 573 - klebsiella pneumoniae'''
SeqProcessor.interact_and_update_dict(filtered_organism_taxid_map)


# In[ ]:


'''Create species_info dictionary that contains keys as species/subspecies/strains and values as sets of taxids'''
species_info = SeqProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)
species_info


# In[ ]:


'''Once again cleaning the names in case we added strains and we dont want to miss valid sequences because of that'''
species_info = {SeqProcessor.filter_key_parts(key): value for key, value in species_info.items()}


# In[ ]:


'''Construct a dictionary from analyzed files that includes seqID as key and taxids as sets of values
add NaN file for entries that had no valid taxids associated with them it's a sign that either blast database or taxonomy db needs updating'''
modified_blast_path = os.path.join(blast_output_path, "modified_blast")
nan_file = os.path.join(modified_blast_path,"NaN acessions.csv")
# Construct dictionary of all taxid - acessions pairs in our data
results_dict = SeqProcessor.construct_dict_from_files(modified_blast_path,nan_file)
'''Most of the time the given sequence is still preserved in the results its presence in the NaN file simply means that one of the acessions got deleted/updated'''


# In[ ]:


'''We can add exceptions now didn't try it with files though probably 
should try to introduce the same ability for interact_and_update_dict'''
exceptions = SeqProcessor.ask_for_exceptions()
print(f"Exceptions to filtering added: {exceptions}")


# In[ ]:


'''Convert values to int'''
for key, value_set in results_dict.items():
    results_dict[key] = {int(val) for val in value_set}    
    
for key, value_set in species_info.items():
    species_info[key] = {int(val) for val in value_set}    


# In[ ]:


exceptions


# In[ ]:


'''This filtering step does the following:
1. Looks through all key-value set pairs in results_dict
2. If one of the values in the set is in exception removes it from the set(but remembers it)
3. If leftover values match atleast one present in species_info 
AND it's the only one associated with a given sequence the record is kept
4. That means that the records which had more than one value associated with it but the other ones were exceptions
are preserved
5.Unfortunetely this didnt help with the Klebsiella_pneumoniae_subsp._pneumoniae completely disappearing from the dataset
6. Still this is not a bug but expected behaviour as KP sequences match to A LOT of taxids'''
filtered_keys = SeqProcessor.filter_with_exceptions(results_dict,species_info,exceptions)

filtered_keys


# In[ ]:


'''Should probably get rid of that but i want to see what happens if we have multi value sets. Did not happen so far.
If errors are encountered at this step just comment out the code and include filtered_keys_final = filtered_keys'''
filtered_keys_final = SeqProcessor.unpack_single_element_sets(filtered_keys)
filtered_keys_final


# In[ ]:


'''Write out to a file a list of rows containing the keys listed in filtered_keys_final'''
SeqProcessor.process_files_with_filter(modified_blast_path,filtered_keys_final)


# In[ ]:


'''Summary function is used here it stores various info about the pipeline - still WIP'''
summary_dir = os.path.join(results_directory,"summary")
SeqProcessor.ensure_directory_exists(summary_dir)
summary_path = os.path.join(summary_dir,"summary.txt")
SeqProcessor.write_summary(summary_path,results_directory)


# In[ ]:


'''Kind of unecessary - there is a simpler way but i have a stable pipeline now and Ill use it to run some experiments'''
fasta_processed_path = SeqProcessor.find_file_by_name('unique_taxids_filtered_sequences_combined_unique_blastn_out.txt',auto= True, folder= results_directory)
fasta_processed_path = Path(fasta_processed_path)
fasta_processed_path = fasta_processed_path.parent
fasta_processed_path


# In[ ]:


'''My naming convention is getting a little bit confusing... if someone wants to fix it go ahead
This function finds twins in our files (the ones with the same number at the end) 
sequences from the same initial sequence keeps them if it find 2'''
SeqProcessor.filter_and_overwrite_files_in_directory(fasta_processed_path,file_pattern="*.txt")


# In[ ]:


unique_directory = SeqProcessor.find_directory_by_name_new('unique_sequences', auto = True, folder=results_directory)


# In[ ]:


sequence_ids_dict, special_dict = SeqProcessor.read_sequence_ids(unique_directory)


# In[ ]:


'''Gets sequence of given id from initial .fasta file (still dealing with short sequences) 
a new function could be written that afterwards strips the numbers from the 'twins' file and uses
the *_results.csv dataframe this way the full sequences can be extracted'''
output_directory = os.path.join(results_directory,"final_output")
filtered_fasta_file = SeqProcessor.find_file_by_name('filtered_sequences_combined_unique.fasta', auto= True, folder= results_directory)
SeqProcessor.ensure_directory_exists(output_directory)
SeqProcessor.filter_fasta_file_dict(filtered_fasta_file,sequence_ids_dict,special_dict,output_directory)


# In[ ]:




