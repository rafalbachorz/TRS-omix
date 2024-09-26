import cffi
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Entrez
import time
import argparse
import sys
from src.SequenceProcessor import SequenceProcessor
from src.FileHandler import FileHandler
from src.pytrsomix import SeqAnalyzer,TRScalculator
from src.BlastProcessor import BLASTProcessor
from src.stats import Stats
import subprocess

#comment

import pathlib
def check_mode(value):
    ivalue = int(value)
    if ivalue < 0 or ivalue > 1:
        raise argparse.ArgumentTypeError(f"Mode must be 0 or 1, got {value}")
    return ivalue

def check_threshold(value):
    fvalue = float(value)
    if fvalue < 0.8 or fvalue > 1.0:
        raise argparse.ArgumentTypeError(f"Identity threshold for cdhit must be between 0.8 and 1.0, got {value}")
    return fvalue

def find_latest_results_directory(base_dir, base_name):
    """
    Find the latest results directory that matches the base name pattern.
    """
    matching_dirs = [d for d in os.listdir(base_dir) if d.startswith(base_name)]
    if not matching_dirs:
        return None
    latest_dir = max(matching_dirs, key=lambda d: os.path.getmtime(os.path.join(base_dir, d)))
    return os.path.join(base_dir, latest_dir)

def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='''This program extracts TRS sequences from a series of input genomes,
                                     allows for length selection of extracted fragments, and prepares sequences for further analysis.''')
    parser.add_argument('--input_fasta_folder_path', help="Path to a folder containing genomes in fasta format from which TRS sequences will be extracted", 
                        required=True)
    parser.add_argument('--tmin', help="Minimum length of TRS sequences", required=True, type=int)
    parser.add_argument('--tmax', help="Maximum length of TRS sequences", required=True, type=int)
    parser.add_argument('--mode', help="Mode of operation, must be 0 or 1", required=True, type=check_mode)
    parser.add_argument('--redo', help="Redo the analysis if results directory already exists", required=False, action='store_true')
    parser.add_argument('--cont', help="Continue the analysis from saved TRS results file", required=False, action='store_true')
    parser.add_argument('--email', help="Address e-mail to be used for connection with NCBI databases", required=True, type=str)
    parser.add_argument('--threshold', help="Identity threshold for clustering using cdhit has to be between 0.8 and 1.0", required=True, type=check_threshold)
    parser.add_argument('--length_to_extract', help="Length of flanking sequences to be extracted from the full TRS sequence", required=True, type=int)
    parser.add_argument('--cd_hit_path', help="Path to the cd-hit-est executable", required=False, type=str)
    parser.add_argument('--taxonomy_db',help='Path to the taxonomy - accession database', required=True, type= str)
    add_ids = parser.add_mutually_exclusive_group()
    add_ids.add_argument('--ids_to_add_to_dictionary', help= "Comma separated list or single value of NCBI IDs", required= False, type=str)
    add_ids.add_argument('--ids_file', help="Path to the file containing NCBI IDs", required= False, type= str)
    exceptions = parser.add_mutually_exclusive_group()
    exceptions.add_argument('--taxids_to_add_to_exceptions', help="Comma separated list or single taxid to add to filtering exceptions", required= False,type=str)
    exceptions.add_argument('--taxids_file', help= "Path to the file containing the taxids that should be added to exceptions", required= False, type=str)
    parser.add_argument('--blast_db', help= "Path to your blast_db", required=True, type=str)
    args = parser.parse_args()

    taxonomy_db = args.taxonomy_db
    ids_to_add_to_dictionary = args.ids_to_add_to_dictionary
    ids_file = args.ids_file
    taxids_to_add_to_exceptions = args.taxids_to_add_to_exceptions
    taxids_file = args.taxids_file



    Entrez.email = SequenceProcessor.validate_and_set_email(args.email)
    print(f"Email adress is currently set to {Entrez.email}")

    # Define the base results directory name and ensure it exists
    input_fasta_folder_path_name = os.path.basename(args.input_fasta_folder_path)
    base_results_directory = os.path.join(os.getcwd(), f"{input_fasta_folder_path_name}_results")
    results_directory = base_results_directory
    FileHandler.ensure_directory_exists(base_results_directory)

    # Define the name of the CSV file that will store the results of the analysis 
    name_of_csv_file_storing_TRS_analysis_results = input_fasta_folder_path_name + "_results.csv"
    path_of_folder_storing_TRS_analysis_results = os.path.join(base_results_directory, "TRS_output")
    FileHandler.ensure_directory_exists(path_of_folder_storing_TRS_analysis_results)
    path_of_csv_file_storing_TRS_analysis_results = os.path.join(path_of_folder_storing_TRS_analysis_results, 
                                                                 name_of_csv_file_storing_TRS_analysis_results)

    if args.cont:
        # Check for the initial results file
        if os.path.exists(path_of_csv_file_storing_TRS_analysis_results):
            combined_trs_results = pd.read_csv(path_of_csv_file_storing_TRS_analysis_results) # no idea how to do that correctly
            print("DataFrame loaded from initial results. Resuming analysis...")
        else:
            # Check for the latest renamed results directory
            latest_results_directory = find_latest_results_directory(os.getcwd(), input_fasta_folder_path_name)
            if latest_results_directory:
                path_of_csv_file_storing_TRS_analysis_results = os.path.join(latest_results_directory, "TRS_output", name_of_csv_file_storing_TRS_analysis_results)
                if os.path.exists(path_of_csv_file_storing_TRS_analysis_results):
                    combined_trs_results = pd.read_csv(path_of_csv_file_storing_TRS_analysis_results)
                    print(f"DataFrame loaded from {latest_results_directory}. Resuming analysis...")
                    results_directory = latest_results_directory
                else:
                    print("No existing results file found to continue from. Exiting program.")
                    sys.exit(1)
            else:
                print("No existing results directory found to continue from. Exiting program.")
                sys.exit(1)
    else:
        # Check if the provided input directory exists
        if not os.path.exists(args.input_fasta_folder_path):
            print("Provided directory path does not exist. Quitting...")
            sys.exit(1)

        print("Existing directory provided. Proceeding....")
        # Check for fasta files in the directory
        fasta_files = [f for f in os.listdir(args.input_fasta_folder_path) if f.endswith('.fasta') or f.endswith('.fa')]
        if not fasta_files:
            print("No fasta files found in the provided directory. Quitting...")
            sys.exit(1)

        # Check if results already exist and handle redo flag
        if os.path.exists(path_of_csv_file_storing_TRS_analysis_results) and not args.redo:
            print("Results already exist. Use --redo to overwrite or --cont to continue from existing results. Exiting program.")
            sys.exit(1)

        print("Starting analysis...")

        trs_calculators = []

        # Iterate over each fasta file and calculate TRS
        for fasta_file in fasta_files:
            path_to_input_fasta = os.path.join(args.input_fasta_folder_path, fasta_file)
            if not os.path.exists(path_to_input_fasta):
                print(f"File '{fasta_file}' does not exist! Skipping....")
                continue

            # Define the TRS file path dynamically (this example assumes a static path, modify as needed)
            trs_file = os.path.join(args.input_fasta_folder_path, 'trs.txt').encode()

            try:
                # Initialize and calculate TRS
                trs_calculator = TRScalculator(sequence=path_to_input_fasta.encode(), trs=trs_file, tmin=args.tmin, tmax=args.tmax, mode=args.mode)
                trs_calculator.calculate()
                trs_calculators.append(trs_calculator)
            except Exception as e:
                print(f"An error occurred while processing '{fasta_file}': {e}")
                continue

        list_of_trs_results = []

        # Iterate over each TRScalculator instance
        for trs_calculator in trs_calculators:
            # Extract results from the calculator
            result = trs_calculator.Result
            # Append the result to the list
            list_of_trs_results.append(result)

        # Concatenate all results into a single DataFrame
        combined_trs_results = pd.concat(list_of_trs_results, ignore_index=True)

        # Remove ">" from >SEQ column
        combined_trs_results[">SEQ"] = combined_trs_results[">SEQ"].str.replace(">","")

        # Save the results of the first analysis step to the CSV file
        combined_trs_results.to_csv(path_of_csv_file_storing_TRS_analysis_results, index=False)
        print(f"Results saved to {path_of_csv_file_storing_TRS_analysis_results}")

    trs_time = time.time()
    print(f"TRS took {trs_time - start_time} seconds")

    l_chars = SequenceProcessor.adjust_input_to_range(args.length_to_extract, args.tmin)
    r_chars = SequenceProcessor.adjust_input_to_range(args.length_to_extract, args.tmin)
    l_chars = int(l_chars)
    r_chars = int(r_chars)
    print(combined_trs_results)
    combined_trs_results = SequenceProcessor.extract_sequences(combined_trs_results, l_chars, r_chars)

    results_directory_after_flanks_extracted = f"{results_directory}_L{l_chars}_R{r_chars}"
    results_directory_after_flanks_extracted_path = os.path.join(os.path.dirname(results_directory), results_directory_after_flanks_extracted)

    print(f"Results directory is set to: {results_directory_after_flanks_extracted_path}")

    if not args.cont:
        if not os.path.exists(results_directory_after_flanks_extracted_path):
            # Rename the existing results directory
            os.rename(results_directory, results_directory_after_flanks_extracted_path)
            print(f"The results directory has been renamed to: {results_directory_after_flanks_extracted_path}")
        else:
            print(f"Directory {results_directory_after_flanks_extracted_path} already exists. Consider using a different name or removing the existing directory.")

    results_directory = results_directory_after_flanks_extracted_path

    ncbi_ids = combined_trs_results["GENOME"].unique().tolist()
    if Entrez.email and Entrez.email == args.email :
         print(f"Email adress is still set to {Entrez.email}")
    organism_map = SequenceProcessor.fetch_organism_names(ncbi_ids,email = Entrez.email)

    # Map NCBI IDs to taxonomic names
    combined_trs_results['Taxonomic Name'] = None
    combined_trs_results['Taxonomic Name'] = combined_trs_results['GENOME'].map(organism_map)
    
    # Identify unmatched genomes
    unmatched_genomes = combined_trs_results[combined_trs_results['Taxonomic Name'].isnull()]["GENOME"].unique()

    if len(unmatched_genomes) > 0:
        print(f"Warning: Some genome IDs could not be matched with taxonomic names: {unmatched_genomes}")
    
    # Extract sequences and create sequence IDs
    combined_trs_results['L_id'] = combined_trs_results['Taxonomic Name'] + '_L' + combined_trs_results['L-No'].astype(str)
    combined_trs_results['R_id'] = combined_trs_results['Taxonomic Name'] + '_R' + combined_trs_results['R-No'].astype(str)
    sequences_df = combined_trs_results[['SEQ_L', 'SEQ_R', 'L_id', 'R_id']]
    print(sequences_df)
    
    path_of_folder_storing_TRS_analysis_results = os.path.join(results_directory, "TRS_output")
    fasta_files_with_flanks = os.path.join(path_of_folder_storing_TRS_analysis_results, "combined_sequences.fasta")
    with open(fasta_files_with_flanks, 'w') as fasta_file:
        for _, row in sequences_df.iterrows():
            # Write left sequence
            fasta_file.write(f'>{row["L_id"]}\n')
            fasta_file.write(f'{row["SEQ_L"]}\n')
            # Write right sequence
            fasta_file.write(f'>{row["R_id"]}\n')
            fasta_file.write(f'{row["SEQ_R"]}\n')

    fasta_files_with_flanks_unique = os.path.join(path_of_folder_storing_TRS_analysis_results, "combined_sequences_unique.fasta")
    SequenceProcessor.rename_sequences(fasta_files_with_flanks, fasta_files_with_flanks_unique)

    cd_hit_results_folder = os.path.join(results_directory, "cd-hit-results")
    FileHandler.ensure_directory_exists(cd_hit_results_folder)
    cd_hit_output_file = os.path.join(cd_hit_results_folder, "combined_sequences_unique_cdhit")
    cd_hit_path = args.cd_hit_path

    results_directory = SequenceProcessor.run_cdhit(cd_hit_path, input_file=fasta_files_with_flanks_unique, output_file=cd_hit_output_file,
                                                    results_directory=results_directory, sc=1, c=args.threshold)
    
    cd_hit_results_folder = os.path.join(results_directory, "cd-hit-results")
    cdhit_clusters_file = os.path.join(cd_hit_results_folder, "combined_sequences_unique_cdhit.clstr")
    non_unique_sequences = os.path.join(cd_hit_results_folder, "combined_sequences_clusters.txt")
    SequenceProcessor.extract_sequence_ids(cdhit_clusters_file, non_unique_sequences)
    clusters_to_be_cleaned = non_unique_sequences
    SequenceProcessor.clean_sequence_ids(clusters_to_be_cleaned)
    clusters_cleaned = clusters_to_be_cleaned
    print(f'current results dir : {results_directory}')
    print(f'{cd_hit_results_folder}')
    fasta_ids_to_remove_because_they_were_in_clusters = FileHandler.read_fasta_ids(clusters_cleaned)
    sequences_after_clusters_filtering_folder = os.path.join(results_directory, "filtered_sequences")
    FileHandler.ensure_directory_exists(sequences_after_clusters_filtering_folder)
    #THERE IS A PROBLEM HERE WITH ASSINGNING VALID FOLDER NAMES

    fasta_with_clustered_ids_removed = os.path.join(sequences_after_clusters_filtering_folder, "not_in_clusters_combined_sequences_unique.fasta")
    fasta_with_clustered_ids_included = os.path.join(sequences_after_clusters_filtering_folder, "in_clusters_combined_sequences_unique.fasta")
    fasta_files_with_flanks_unique =os.path.join(results_directory,"TRS_output","combined_sequences_unique.fasta")
    blast_results_folder = os.path.join(results_directory, "blast_output")
    FileHandler.ensure_directory_exists(blast_results_folder)
    FileHandler.filter_fasta_file(fasta_files_with_flanks_unique, fasta_with_clustered_ids_removed, fasta_ids_to_remove_because_they_were_in_clusters)
    FileHandler.filter_fasta_file_clusters(fasta_files_with_flanks_unique, fasta_with_clustered_ids_included, fasta_ids_to_remove_because_they_were_in_clusters)
    
    # Paths to your SLURM and fallback script
    slurm_script = 'TRS_BLAST.sh'
    fallback_script = 'TRS_BLAST_NOSLURM.sh'

    # Set the paths to the query directory and BLAST database
    query_dir = sequences_after_clusters_filtering_folder  # Replace with your actual path
    blast_db = args.blast_db # Replace with your actual path

    # Set the environment variables before running the SLURM script
    env = os.environ.copy()
    env['QUERY_DIR'] = query_dir
    env['BLAST_DB'] = blast_db
    env['OUTPUT_DIR'] = blast_results_folder

    try:
        # Attempt to submit the SLURM job using sbatch with the --wait option
        subprocess.run(['sbatch', '--wait', slurm_script], env=env, check=True)

        # If the SLURM job completes successfully, continue with the Python script
        print("SLURM job completed successfully. Continuing with the Python script.")

    except subprocess.CalledProcessError as e:
        print(f"SLURM job submission failed with error: {e}. Running fallback script...")

        # If SLURM job submission fails, run the fallback script and wait for its completion
        result = subprocess.run([fallback_script, query_dir, blast_db], check=True)

        if result.returncode == 0:
            print("Fallback script completed successfully. Continuing with the Python script.")
        else:
            print("Fallback script failed.", file=sys.stderr)
            sys.exit(1)
    
    #BLAST part
    blast_output_path = blast_results_folder
    FileHandler.convert_to_txt(blast_output_path)
    FileHandler.filter_and_overwrite_blast_file(blast_output_path)
    modified_blast_path = os.path.join(blast_output_path, "modified_blast")
    FileHandler.ensure_directory_exists(modified_blast_path)
    

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

    Entrez.email = SequenceProcessor.validate_and_set_email(args.email)
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
