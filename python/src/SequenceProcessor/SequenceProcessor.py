import os
import re
import subprocess
import pandas as pd
from Bio import Entrez, SeqIO
from tqdm import tqdm
import time


class SequenceProcessor:
    @staticmethod
    def calculate_sequence_lengths(tmin):
        """
        Calculate maximum lengths for left and right sequences based on a minimum total length.
        """
        try:
            l_chars_max = tmin // 2
            r_chars_max = tmin // 2
            print(f"The maximum available length for SEQ_L and SEQ_R is {l_chars_max} and {r_chars_max} respectively")
            return l_chars_max, r_chars_max
        except Exception as e:
            print(f"An error has occurred: {e}")
    @staticmethod
    def adjust_input_to_range(length_to_extract,tmin): # to be rewriten
        """
        Prompts the user for input and adjusts it to ensure it does not exceed the maximum allowed value.
        """
        while True:
            try:
                max_val = tmin/2
                length_to_extract = length_to_extract
                if length_to_extract > max_val:
                    print(f"Your input was adjusted to {max_val}.")
                    return max_val
                return length_to_extract
            except ValueError:
                print("Please ensure that you provide a valid integer")
            except Exception as e:
                print(f"An unexpected error has occurred: {e}")
    @staticmethod
    def extract_sequences(combined_results, l_chars, r_chars):
        """
        Extracts left and right sequences from the combined results based on the specified lengths.
        """
        try:
            combined_results['SEQ_L'] = combined_results['>SEQ'].str.slice(0,l_chars)
            combined_results['SEQ_R'] = combined_results['>SEQ'].str.slice(-r_chars)
            print("Dataframe filtered successfully.")
            return combined_results
        except Exception as e:
            print(f"An error has occured: {e}!")
    @staticmethod
    def validate_and_set_email(email):
        """
        Validates inputted email address, sets it for Entrez access.
        """
        def validate_email(email):
            return bool(re.match(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$", email))
        try:
            if validate_email(email):
                            Entrez.email = email
                            print("Email address set successfully")
            else:
                print("Invalid email format. Please enter a valid email address.")
        except Exception as e:
             print(f"An unexpected error has occurred {e}. Please try again")
        return Entrez.email
             
    @staticmethod
    def fetch_organism_names(ncbi_ids, max_retries = 3, rate_limit = 0.33,email=Entrez.email,timeout =30):
        """
        Fetch organism names from NCBI based on NCBI IDs.
        """
        Entrez.email = email
        organism_map = {}
        for ncbi_id in ncbi_ids:
            start_time = time.time()
            retries = 0
            success = False
            while retries < max_retries and not success:
                try:
                    print(f"Fetching organism name for id: {ncbi_id}")
                    handle = Entrez.efetch(db='nucleotide',id=ncbi_id,rettype ='gb',retmode = 'text')
                    record = handle.read()
                    for line in record.splitlines():
                        if line.startswith("  ORGANISM"): #Formating quirk
                                organism = line.split("ORGANISM")[1].strip().replace(" ","_")
                                organism_map[ncbi_id] = organism
                                if organism != 'nan':
                                    print(f"Retrieved organism name for id: {ncbi_id} successfully {organism}.")
                                    success = True
                                else:
                                    raise Exception('Organism names is nan')
                except Exception as e:
                    print(f"Error fetching organism name for NCBI ID: {ncbi_id}: {e}! Retrying...")
                    retries += 1
                    time.sleep(rate_limit)
                finally:
                    if not success and retries >= max_retries:
                        time.sleep(rate_limit)
            
            elapsed_time = time.time() - start_time
            
            if elapsed_time > timeout:
                print(f"Timeout exceeded. Restarting from ncbi_id {ncbi_id}.")
                return SequenceProcessor.fetch_organism_names(ncbi_ids,email=Entrez.email)
        return organism_map
            
    @staticmethod
    def rename_sequences(input_file,output_file):
        """
        Renames sequences in a FASTA file.
        """
        try:
            with open(input_file,"r") as input_handle, open(output_file,"w") as output_handle:
                    records = list(SeqIO.parse(input_handle,"fasta"))
                    input_handle.seek(0)
                    pair_index = 1
                    for record in tqdm(SeqIO.parse(input_handle,"fasta"), total=len(records), desc="Renaming sequences"):
                        original_id = record.id
                        new_id = f"{original_id}_{pair_index}"
                        record.id = new_id
                        record.description = ""
                        SeqIO.write(record,output_handle,"fasta")
                        if 'R' in original_id:
                                pair_index += 1
                    print(f"Sequences successfully renamed and saved at {output_file}!")
        except FileNotFoundError:
            print(f"Error: The file {input_file} was not found!")
        except PermissionError:
            print(f"Permission denied when accessing the files!")
        except IOError as e:
            print(f"An I/O error occurred: {e}")
        except Exception as e:
            print(f"An unexpected error has occurred: {e}")
    @staticmethod
    def adjust_word_length(c):
        """
        Dynamically adjusts 'n' based on sequence identity threshold 'c'
        """
        try:
            threshold_n_values = [
                (0.95,10),
                (0.9,9),
                (0.88,8),
                (0.85,7),
                (0.8,6),
                (0.75,5)
            ]
            if c==1.0:
                return 11
            for lower_bound, n_value in threshold_n_values:
                if c>= lower_bound:
                    print(f"Wordsize was set to {n_value} for {c}")
                    return n_value
            return None
        except TypeError as e:
            print(f"Error: {e}")
        except ValueError as e:
            print(f"Error: {e}")
        except Exception as e:
            print(f"An unexpected error has occurred: {e}")
    @staticmethod
    def run_cdhit(
         cdhit_path,
         input_file,
         output_file,
         c=None,
         d=0,
         m="t",
         g=0,
         G=1,
         sc=None,
         results_directory="results" #wtf ??
    ):
        """
        Run cdhit directly from python adjusting n based on the user provided c
        """
        n = SequenceProcessor.adjust_word_length(c)

        if not os.path.exists(os.path.join(cdhit_path, "cd-hit-est")):
                raise FileNotFoundError("CD-HIT executable not found in the specified directory.")
        from ..FileHandler.FileHandler import FileHandler
        FileHandler.ensure_directory_exists(results_directory)

        cmd = [
             os.path.join(cdhit_path,'cd-hit-est'),
             "-i", input_file,
             "-o", output_file,
             "-c", str(c),
             "-n", str(n),
             "-d", str(d),
             "-M", m,
             "-g", str(g),
             "-G", str(G),
        ]
        
        if sc is not None:
             cmd.extend(["-sc", str(sc)])
        try:
            subprocess.run(cmd, check=True)
            new_results_directory = f"{results_directory}_c{c}"
            if not os.path.exists(new_results_directory):
                os.rename(results_directory, new_results_directory)
                results_directory = new_results_directory
                print(f"The results directory has been renamed to: {results_directory}")
            else:
                print(f"Warning: The directory {new_results_directory} already exists. Results directory was not renamed.")
                results_directory = new_results_directory
        except subprocess.CalledProcessError as e:
            print(f"CD-HIT command failed: {e}")
            raise
        
        return results_directory
    
    @staticmethod
    def filter_taxonomy_file(taxonomy_file,accessions,chunksize= 50000):
        """
        Filters the taxonomy file in chunks, keeping the rows matching accessions
        """
        try:
            filtered_rows = []
            cols_to_use = [1,2]
            chunk_count = 0
            estimated_chunks = 6548

            for chunk in tqdm(pd.read_csv(taxonomy_file,sep="\t",usecols=cols_to_use, chunksize= chunksize), total= estimated_chunks):
                filtered_chunk = chunk[chunk.iloc[:,0].isin(accessions)]   
                filtered_rows.append(filtered_chunk)
            filtered_df = pd.concat(filtered_rows,ignore_index=True)

            return filtered_df
        
        except FileNotFoundError:
            print(f"Error: The file {taxonomy_file} was not found.")
        except PermissionError:
            print(f"Error access to the file {taxonomy_file} was denied.")
        except pd.errors.EmptyDataError:
            print(f"Error: The file {taxonomy_file} is empty or all rows have been filtered out")
        except Exception as e:
            print(f"An unexpected error has occurred: {e}")
    
    @staticmethod
    def process_sequences(input_file, output_file, processor_func):
        """
        General function to process sequences using a provided processor function.
        """
        try:
            with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    processed_record = processor_func(record)
                    SeqIO.write(processed_record, output_handle, "fasta")
        except FileNotFoundError:
            print(f"Error: The file {input_file} was not found!")
        except PermissionError:
            print(f"Permission denied when accessing the files!")
        except IOError as e:
            print(f"An I/O error occurred: {e}")
        except Exception as e:
            print(f"An unexpected error has occurred: {e}")

    @staticmethod
    def rename_sequences(input_file, output_file):
        """
        Renames sequences in a FASTA file using a general processor function.
        """
        def rename_processor(record):
            nonlocal pair_index
            original_id = record.id
            new_id = f"{original_id}_{pair_index}"
            record.id = new_id
            record.description = ""
            if 'R' in original_id:
                pair_index += 1
            return record

        pair_index = 1
        SequenceProcessor.process_sequences(input_file, output_file, rename_processor)
    
    @staticmethod
    def fetch_organism_taxids(ncbi_ids,rate_limit = 0.33, max_retries = 3, timeout = 30):
        """
        Fetch organism TaxIDs from NCBI based on NCBI IDs, using organism names as a fallback.
        """
        organism_taxid_map = {}
        organism_names = SequenceProcessor.fetch_organism_names(ncbi_ids)

        for ncbi_id in ncbi_ids:
            retries = 0
            success = False
            while retries < max_retries and not success:
                try:
                    handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
                    record = handle.read()
                    handle.close()
                    taxid = None
                    for line in record.splitlines():
                        if line.strip().startswith("/db_xref=\"taxon:"): 
                            taxid = line.split(":")[1].replace("\"", "").strip()
                            break
                    if taxid != 'NaN':
                        organism_taxid_map[ncbi_id] = taxid
                        print(f"Retrieved TaxID {taxid} for NCBI ID {ncbi_id} successfully.")
                        success = True
                    else:
                        print(f"TaxID not found for {ncbi_id}")
                except Exception as e:
                    print(f"Error fetching TaxID for NCBI ID {ncbi_id}: {e}. Retrying....")
                    retries += 1
                    time.sleep(rate_limit)
                finally:
                    if not success and retries >= max_retries:
                        time.sleep(rate_limit)

            organism_taxid_map_with_names = {
                organism_names.get(ncbi_id, "Unknown_Organism"): taxid
                for ncbi_id, taxid in organism_taxid_map.items()
            }

        print("All organism TaxIDs successfully retrieved and mapped to organism names.")
        return organism_taxid_map_with_names
    
    @staticmethod #SHOULD BE MOVED TO BLASTPROCESSOR
    def filter_key_parts(key):
        """
        Filters parts of the key based on specific rules.
        """
        parts = key.split("_")
        filtered_parts = [parts[0]]  # Always include the first part

        for part in parts[1:]:
            if not re.search(r'^[A-Z]', part) and not re.search(r'\d', part):
                filtered_parts.append(part)
        
        return '_'.join(filtered_parts)
    
    @staticmethod
    def process_taxid_input(input_str):
        """
        Process the input string which could be single taxid, comma separated list of taxids or taxid file with one taxid per row.
        """
        taxids = set()
        for taxid in input_str.split(','):
            taxid = taxid.strip()
            if taxid.isdigit():
                taxids.add(taxid)
            else:
                print(f"Warning: Badly formatted input detected ('{taxid}').")
        return taxids
    @staticmethod
    def count_fasta_sequences(file_path):
        """
        Counts the number of sequences in a FASTA file.
        """
        try:
            with open(file_path, "r") as fasta_file:
                return sum(1 for line in fasta_file if line.startswith('>'))
        except FileNotFoundError:
            print(f"File {file_path} not found.")
            return None
    
    @staticmethod
    def count_clusters(file_path):
        """
        Counts the total number of clusters and the number of clusters with more than 2 entries in a CD-HIT cluster file.
        """
        total_clusters = 0
        clusters_with_more_than_two_entries = 0
        current_cluster_size = 0

        try:
            with open(file_path, "r") as file:
                for line in file:
                    if line.startswith('>Cluster'):
                        if current_cluster_size > 2:
                            clusters_with_more_than_two_entries += 1
                        total_clusters += 1
                        current_cluster_size = 0
                    else:
                        current_cluster_size += 1

                if current_cluster_size > 2:
                    clusters_with_more_than_two_entries += 1
        
        except FileNotFoundError:
            print(f"File {file_path} not found")
            return None, None
        
        return total_clusters, clusters_with_more_than_two_entries
    
    @staticmethod
    def count_lines_in_file(file_path):
        """
        Counts the number of lines in a file.
        """
        try:
            with open(file_path, "r") as file:
                return sum(1 for _ in file)
        except FileNotFoundError:
            print(f"File {file_path} not found.")
            return None
        
    @staticmethod
    def extract_sequence_ids(input_file, output_file):
        """
        Extract sequence IDs from clusters with more than 2 sequences and save them to a file.

        Parameters:
            input_file (str): Path to the input cluster data file.
            output_file (str): Path to the output file to save sequence IDs.
        """
        try:
            with open(input_file, "r") as f, open(output_file, "w") as outfile:
                current_cluster = []
                for line in f:
                    if line.startswith(">Cluster"):
                        if len(current_cluster) >= 2:
                            for seq in current_cluster:
                                outfile.write(seq)
                        current_cluster = []  # Reset for next cluster
                    elif line.strip():  # Non-empty line
                        current_cluster.append(line)  # Add sequence ID to current cluster
                # Check last cluster in file
                if len(current_cluster) >= 2:
                    for seq in current_cluster:
                        outfile.write(seq)
        except FileNotFoundError:
            print(f"Error: The file {input_file} does not exist.")
        except PermissionError:
            print(f"Error: Permission denied when acessing {input_file} or writing to {output_file}!")
        except Exception as e:
            print(f"An unexpected error has occured.")
    @staticmethod
    def clean_sequence_ids(input_file):
        """
        Clean sequence IDs file by removing prefixes and suffixes.

        Parameters:
            input_file (str): Path to the input file with sequence IDs.

        Returns:
            None
        """
        try:
            with open(input_file, "r") as infile:
                lines = infile.readlines()

            with open(input_file, "w") as outfile:
                for line in tqdm(lines, desc="Cleaning sequence IDs"):
                    sequence_id = line.split(">")[1].split("...")[0].strip()
                    outfile.write(">" + sequence_id + "\n")
        except FileNotFoundError:
            print(f"Error: file {input_file} does not exist!")
        except PermissionError:
            print(f"Error: Access to {input_file} denied")
        except Exception as e:
            print(f"An unexpected error has occurred: {e}.")

        print(f"Overwriting file {input_file} with cleaned entries.")
