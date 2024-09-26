# %%
import cffi
import os
import numpy as np
import pandas as pd
from io import StringIO
from enum import Enum
from Bio import SeqIO
from Bio import Entrez
from pathlib import Path
import matplotlib.pyplot as plt
from multiprocessing import Pool
import glob
import re
import fnmatch
from src.SequenceProcessor import SequenceProcessor
from src.FileHandler import FileHandler
from src.pytrsomix import SeqAnalyzer,TRScalculator
from src.BlastProcessor import BLASTProcessor
import argparse
import time

# %%
def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='''This program analyzes the blast output files to find species specific sequences coming from genomes
                                     analyzed in the previous step''')
    parser.add_argument('--blast_output_path', help="Path to a folder containing blast results in a valid format", type=str, required=True)
    parser.add_argument('--taxonomy_db',help='Path to the taxonomy - accession database', required=True, type= str)
    parser.add_argument('--email', help= "Addres email to be used for connection with NCBI servers", required= True, type= str)
    add_ids = parser.add_mutually_exclusive_group()
    add_ids.add_argument('--ids_to_add_to_dictionary', help= "Comma separated list or single value of NCBI IDs", required= False, type=str)
    add_ids.add_argument('--ids_file', help="Path to the file containing NCBI IDs", required= False, type= str)
    exceptions = parser.add_mutually_exclusive_group()
    exceptions.add_argument('--taxids_to_add_to_exceptions', help="Comma separated list or single taxid to add to filtering exceptions", required= False,type=str)
    exceptions.add_argument('--taxids_file', help= "Path to the file containing the taxids that should be added to exceptions", required= False, type=str)
    args = parser.parse_args()

    blast_output = args.blast_output_path
    taxonomy_db = args.taxonomy_db
    email = args.email
    ids_to_add_to_dictionary = args.ids_to_add_to_dictionary
    ids_file = args.ids_file
    taxids_to_add_to_exceptions = args.taxids_to_add_to_exceptions
    taxids_file = args.taxids_file

    #FileHandler.convert_to_txt(blast_output)
    
    #FileHandler.filter_and_overwrite_blast_file(blast_output)
    blast_output_path = Path(blast_output)
    results_directory = blast_output_path.parent
    modified_blast_path = os.path.join(blast_output_path, "modified_blast")
    FileHandler.ensure_directory_exists(modified_blast_path)

    #tax_df = pd.read_csv(taxonomy_db,sep="\t",usecols = [1,2])
    #print("Initial Taxonomy DataFrame:\n", tax_df.head(10))
    #Collect accessions from BLAST files
    accessions = BLASTProcessor.collect_accessions_from_blast_files(blast_output_path)
    
    # Filter taxonomy file based on collected accessions
    tax_df = SequenceProcessor.filter_taxonomy_file(taxonomy_db, accessions)
    print("Collected Accessions:", accessions)
    print("Filtered Taxonomy DataFrame:\n", tax_df)

    # Create a mapping from accession to taxid
    print(f"Creating a mapping between accessions and taxids....")
    taxid_accessions_dict = {}
    for index, row in tax_df.iterrows():
        accession = row[tax_df.columns[0]]
        taxid = row[tax_df.columns[1]]

        if taxid in taxid_accessions_dict:
            taxid_accessions_dict[taxid].append(accession)
        else:
            taxid_accessions_dict[taxid] = [accession]

    accession_to_taxid = {accession: taxid for taxid, accessions in taxid_accessions_dict.items() for accession in accessions}

    # Process BLAST files to map accessions to TaxIDs and add them as a new column
    print(f"Appending a new column with taxids to the blast files.....")
    BLASTProcessor.match_accessions_to_taxids_and_add_column(
        blast_output_path,
        modified_blast_path,
        lambda accession: BLASTProcessor.map_accession_to_taxid(accession, accession_to_taxid)
    )

    end_time = time.time()
    print(f"Processing completed in {end_time - start_time:.2f} seconds")
    print(f"Searching for results file...")
    results_file = FileHandler.find_file_by_name('*_results.csv', folder= results_directory)
    if results_file:
        result_file = results_file[0]
        print(f"Results file found at : {results_file}")
    else:
        print("No '*_results.csv' file selected or found.")
    
    combined_results = pd.read_csv(result_file)

    Entrez.email = SequenceProcessor.validate_and_set_email(email)
    print(f"Current email is set to {Entrez.email}")
    ncbi_ids = combined_results["GENOME"].unique().tolist()
    tax_map = SequenceProcessor.fetch_organism_taxids(ncbi_ids)

    filtered_organism_taxid_map = {SequenceProcessor.filter_key_parts(key): value for key, value in tax_map.items()}
    print(f"Species - taxid pairs detected in dataset : {filtered_organism_taxid_map}")

    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)

    BLASTProcessor.interact_and_update_dict(filtered_organism_taxid_map,ncbi_ids_list=ids_to_add_to_dictionary,file_path=ids_file)

    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)
    species_info = {SequenceProcessor.filter_key_parts(key): value for key, value in species_info.items()}

    nan_file = os.path.join(modified_blast_path,"NaN acessions.csv")

    # Construct dictionary of all taxid - acessions pairs in our data
    results_dict = BLASTProcessor.construct_dict_from_files(modified_blast_path,nan_file)

    exceptions = BLASTProcessor.ask_for_exception(exception_ids=taxids_to_add_to_exceptions,file_path=taxids_file)
    if exceptions:
        print(f"Exceptions to filtering added: {exceptions}")
    else:
        print(f"No exceptions provided!")
    
    '''Convert values to int'''
    for key, value_set in results_dict.items():
        results_dict[key] = {int(val) for val in value_set}    
        
    for key, value_set in species_info.items():
        species_info[key] = {int(val) for val in value_set}  

    '''This filtering step does the following:
    1. Looks through all key-value set pairs in results_dict
    2. If one of the values in the set is in exception removes it from the set(but remembers it)
    3. If leftover values match atleast one present in species_info 
    AND it's the only one associated with a given sequence the record is kept
    4. That means that the records which had more than one value associated with it but the other ones were exceptions
    are preserved
    5.Unfortunetely this didnt help with the Klebsiella_pneumoniae_subsp._pneumoniae completely disappearing from the dataset
    6. Still this is not a bug but expected behaviour as KP sequences match to A LOT of taxids'''
    
    print("Filtering keys....")
    filtered_keys = BLASTProcessor.filter_with_exceptions(results_dict,species_info,exceptions)
    filtered_keys_final = BLASTProcessor.unpack_single_element_sets(filtered_keys)

    print("Filtering files...")
    FileHandler.process_files_with_filter(modified_blast_path,filtered_keys_final)
    
    processed_fasta_file_path = FileHandler.find_file_by_name(file_name="unique_taxids_in_clusters_combined_sequences_unique_blastn_out.txt",folder=results_directory)
    print(f"{processed_fasta_file_path}")
    processed_fasta_file_path = processed_fasta_file_path[0]
    processed_fasta_file_path = Path(processed_fasta_file_path)
    processed_fasta_file_path = processed_fasta_file_path.parent
    print(f"{processed_fasta_file_path}")


    BLASTProcessor.separate_into_singles_and_twins(processed_fasta_file_path,file_pattern="*.txt")
    sequence_ids_dict, special_dict = BLASTProcessor.read_sequence_ids(processed_fasta_file_path)


    filtered_fasta_file = FileHandler.find_file_by_name('not_in_clusters_combined_sequences_unique.fasta',folder= results_directory)
    filtered_fasta_file = filtered_fasta_file[0]
    filtered_fasta_file = Path(filtered_fasta_file)
    output_directory = os.path.join(results_directory,"final_output")
    FileHandler.ensure_directory_exists(output_directory)
    fasta_to_parse = BLASTProcessor.filter_fasta_file_dict(filtered_fasta_file,sequence_ids_dict,special_dict,output_directory)
    BLASTProcessor.extract_full_TRS_sequences(combined_results,fasta_to_parse,results_directory)
    


if __name__ == "__main__":
    main()



