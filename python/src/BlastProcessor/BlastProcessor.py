import os
import fnmatch
import pandas as pd
from src.SequenceProcessor import SequenceProcessor
from Bio import Entrez
from pathlib import Path

class BLASTProcessor:

    @staticmethod
    def process_blast_files_in_directory(directory_path, file_pattern="*.txt"):
        """
        Processes BLAST output files in a specified directory.
        
        Args:
            directory_path (str): Path to the directory containing BLAST files.
            file_pattern (str): Pattern to match BLAST files (default is "*.txt").
        """
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                BLASTProcessor.filter_and_overwrite_blast_file(file_path)

    @staticmethod
    def filter_and_overwrite_blast_file(file_path):
        """
        Processes BLAST output file to remove duplicate pairs of sequence IDs and accession numbers,
        keeps the first occurrence of each unique pair along with the entire line of data.
        
        Args:
            file_path (str): Path to the BLAST file to be processed.
        """
        unique_pairs = set()
        lines_to_keep = []

        try:
            with open(file_path, "r") as file:
                for line in file:
                    try:
                        columns = line.strip().split('\t')
                        if len(columns) >= 2:
                            pair = (columns[0], columns[1])
                            if pair not in unique_pairs:
                                unique_pairs.add(pair)
                                lines_to_keep.append(line)
                    except Exception as e:
                        print(f"An unexpected error has occurred while processing a line in {file_path}: {e}")
                print(f"Searching for unique sequence ID - accession pairs in {file_path}....")
        except PermissionError as e:
            print(f"Permission denied while trying to access {file_path}: {e}")
            return
        except IOError as e:
            print(f"An error has occurred while opening or reading {file_path}: {e}")
            return
        
        try:
            with open(file_path, "w") as file:
                file.writelines(lines_to_keep)
        except PermissionError as e:
            print(f"Permission denied while accessing {file_path} for writing: {e}")
        except IOError as e:
            print(f"An error occurred while writing to {file_path}: {e}")

    @staticmethod
    def collect_accessions_from_blast_files(directory_path, file_pattern="*.txt"):
        """
        Collects accession numbers from BLAST files in a specified directory.
        
        Args:
            directory_path (str): Path to the directory containing BLAST files.
            file_pattern (str): Pattern to match BLAST files (default is "*.txt").
            
        Returns:
            set: A set of unique accession numbers.
        """
        accessions = set()
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                with open(file_path, "r") as file:
                    for line in file:
                        columns = line.strip().split('\t')
                        if len(columns) >= 2:
                            accessions.add(columns[1])
        return accessions

    @staticmethod
    def filter_and_overwrite_files_in_directory(directory_path, file_pattern="*.txt"):
        """
        Filters and overwrites files in a directory based on BLAST results.
        
        Args:
            directory_path (str): Path to the directory containing files to filter.
            file_pattern (str): Pattern to match files (default is "*.txt").
        """
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                BLASTProcessor.filter_and_overwrite_blast_file(file_path)
    @staticmethod
    def map_accession_to_taxid(accession, taxid_dict):
        return taxid_dict.get(accession, '')
    @staticmethod
    def match_accessions_to_taxids_and_add_column(blast_output_path, modified_blast_path, map_func):
        for filename in os.listdir(blast_output_path):
            if filename.endswith(".txt"):
                input_file_path = os.path.join(blast_output_path, filename)
                
                # Read the file into a DataFrame
                df = pd.read_csv(input_file_path, sep='\t', header=None)  # No column names specified
                
                # Map accessions to TaxIDs and add them as a new column after the last existing column
                df[df.shape[1]] = df[1].map(lambda x: map_func(x))  # Apply the mapping function to column 1
                
                # Construct output file path
                modified_file_path = os.path.join(modified_blast_path, f"taxids_{filename}")
                
                # Save the modified DataFrame to a new file in the output directory
                df.to_csv(modified_file_path, sep='\t', index=False, header=False)
    @staticmethod
    def append_taxids_to_filtered_map(filtered_map):
        """
        For each organism name in the filtered map, query NCBI to find the corresponding TaxID and add it as
        another value to the key, convert dictionary to use sets for values. If value found is already present skips the ID

        Args:
        - filtered_map (dict) - dictionary mapping filtered organisms to their TaxIDs
        """
        updated_map = {}
        for organism_name, original_taxid in filtered_map.items():
            # Makes is so the values are stored in sets 
            if organism_name not in updated_map:
                updated_map[organism_name] = set()
            updated_map[organism_name].add(original_taxid)

            try:
                handle = Entrez.esearch(db="taxonomy", term=organism_name, retmode="xml")
                search_results = Entrez.read(handle)
                handle.close()
                taxids = search_results.get("IdList",[])
                if taxids:
                    #Assume the first result is the most relevant
                    new_taxid = taxids[0]
                    if new_taxid not in updated_map[organism_name]:
                        updated_map[organism_name].add(new_taxid)
                        print(f"Updated {organism_name} with additional TaxID: {taxids[0]}")
                    else:
                        print(f"TaxID {new_taxid} already present for {organism_name}, skipping....")
                else:
                    print(f"No additional TaxID found for: {organism_name}")
            except Exception as e:
                print(f"Error fetching additional TaxID for {organism_name}: {e}")
        return updated_map
    @staticmethod
    def interact_and_update_dict(dictionary, ncbi_ids_list = None, file_path = None):
        """
        Interacts with user to possibly add more TaxIDs and organism names to species_info dictionary allows to specify a file_path containing 
        NCBI ids(one per line) to include find and include in the dictionary

        Args:
        dictionary(dict): The existing dictionary mapping organism names to SETS of TaxIDs.
        user_input(str): Comma separated list or single value of NCBI IDs.
        file_path(str): Path to file containing NCBI IDs.
        """
        if ncbi_ids_list:
            ids = [id.strip() for id in ncbi_ids_list.split(",")]
            if ids:
                new_organism_ids = SequenceProcessor.fetch_organism_taxids(ids)
                for organism_name, taxid_set in new_organism_ids.items():
                    if organism_name in dictionary:
                        print(f"Organism {organism_name} already present in dictionary,updating dictionary")
                        dictionary[organism_name].update(taxid_set)
                    else:
                        print(f"New organism detected : {organism_name}")
                        dictionary[organism_name] = taxid_set
        elif file_path:
            try:
                with open(file_path,'r') as file:
                    ids = [line.strip() for line in file if line.strip()]
                    new_organism_ids = SequenceProcessor.fetch_organism_taxids(ids)
                    for organism_name, taxid_set in new_organism_ids.items():
                        if organism_name in dictionary:
                            print(f"Organism {organism_name} already present in dictionary,updating dictionary")
                            dictionary[organism_name].update(taxid_set)
                        else:
                            print(f"New organism detected : {organism_name}")
                            dictionary[organism_name] = taxid_set
            except FileNotFoundError:
                print("The specified file was not found, check path and try again")
            except Exception as e:
                print("An unexpected error has occured: {e}")
        else:
            print("No input provided nothing will be added")
    @staticmethod
    def ask_for_exception(exception_ids = None, file_path = None):
        """
        Adds exceptions to further filtering steps based on provided arguments.

        Args: 
            exception_ids(str): Comma separated list or single value of taxids.
            file_path(str): Path to file containing the exceptions.
        """
        exceptions = set()

        if exception_ids:
            exceptions = SequenceProcessor.process_taxid_input(exception_ids)
        elif file_path:
            try:
                with open(file_path, 'r') as file:
                    file_exceptions = set(line.strip() for line in file if line.strip().isdigit())
                    exceptions.update(file_exceptions)
            except FileNotFoundError:
                print("The specified file was not found, check path and try again")
            except Exception as e:
                print(f"An unexpected error has occurred: {e}")
        else:
            print("No input provided. No exceptions will be added.")

        return exceptions
    
    @staticmethod
    def filter_with_exceptions(results_dict, species_info, exceptions):
        """
        Filters the results dictionary based on the presence of taxids in the species info dictionary, taking into account exceptions.
        """
        filtered_keys = {}
        all_species_taxids = set.union(*species_info.values())  # Combine all taxids from species info for easier lookup

        for key, values in results_dict.items():
            non_exception_values = values - exceptions  # Remove exceptions 
            if non_exception_values:  # If there are none left
                if len(non_exception_values) == 1 and non_exception_values & all_species_taxids:
                    filtered_keys[key] = non_exception_values
            else:  # If all values are exceptions 
                if values & all_species_taxids:
                    filtered_keys[key] = values
        return filtered_keys
    @staticmethod
    def unpack_single_element_sets(input_dict):
        """
        Processes the input dictionary to replace any single-element set values with the element itself.

        Args:
        - input_dict (dict): The dictionary whose single-element sets are to be unpacked.

        Returns:
        - dict: A dictionary where each single-element set value has been replaced by the element itself.
        """
        processed_dict = {}

        for key, value in input_dict.items():
            # Check if the value is a set with exactly one entry
            if isinstance(value, set) and len(value) == 1:
                # "Unpack" the set, storing its single element as the value
                processed_dict[key] = next(iter(value))
            else:
                # If the value is not a single-element set, retain it as is
                processed_dict[key] = value

        return processed_dict
    @staticmethod
    def separate_into_singles_and_twins(directory_path, file_pattern="*.txt"):
        """
        Filters each file in a specified directory based on the last number in the identifiers of its lines,
        then overwrites each file with its filtered content for entries with paired identifiers.
        Writes entries with unique identifiers to a separate file named '<original_filename>_singles.txt'.
        Targets files matching a given file pattern (default is "*.txt").

        Args:
            directory_path (str): Path to the directory containing the files to process.
            file_pattern (str): Pattern of the file names to process (default is "*.txt").
        """
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                single_entries_path = os.path.join(directory_path, f"{os.path.splitext(filename)[0]}_singles.txt")

                print(f"Processing file: {filename}")  # Debugging

                with open(file_path, 'r') as file:
                    lines = file.readlines()

                print(f"Total lines read from {filename}: {len(lines)}")  # Debugging

                number_count = {}
                filtered_entries = []
                non_unique_entries = []

                # First pass to count occurrences of each identifier
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        number_count[number] = number_count.get(number, 0) + 1

                print(f"Number count: {number_count}")  # Debugging

                # Second pass to separate unique and non-unique entries
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        if number_count[number] > 1:
                            filtered_entries.append(line.strip())
                        else:
                            non_unique_entries.append(line.strip())

                print(f"Filtered entries for {filename}: {len(filtered_entries)}")  # Debugging
                print(f"Non-unique entries for {filename}: {len(non_unique_entries)}")  # Debugging

                # Overwrite the original file with paired identifiers
                with open(file_path, 'w') as file:
                    for entry in filtered_entries:
                        file.write(entry + '\n')

                print(f"Overwritten {filename} with filtered entries.")  # Debugging

                # Write out single entries to a separate file
                with open(single_entries_path, 'w') as file:
                    for entry in non_unique_entries:
                        file.write(entry + '\n')

                print(f"Written non-unique entries to {single_entries_path}.")  # Debugging

    @staticmethod
    def construct_dict_from_files(directory, nan_file):
        '''
        Construct dictionaries containing the sequence ids as keys and taxids as values

        Args:
        -directory (str) - path to files from which to construct dictionary 
        '''
        data_dict = {}
        nan_keys = set()
        for filename in os.listdir(directory):
            if filename.endswith('.txt'):
                filepath = os.path.join(directory,filename)
                # determine the number of columns
                with open(filepath,"r") as f:
                    first_line = f.readline()
                    num_columns = len(first_line.split('\t'))
                
                #Read only necessary columns (1st and last)
                df = pd.read_csv(filepath, sep='\t',header=None, usecols=[0,num_columns-1])

                #Group by the first column and convert the last one to sets
                grouped = df.groupby(0)[num_columns-1].apply(set).to_dict()

                #Merge the current file's dictionary with the main dictionary
                for key, value_set in grouped.items():
                    '''
                    This is a bandaid fix for a problem that can be encountered when either taxonomy or blast database is out of date 
                    '''
                    value_set = {val for val in value_set if pd.notna(val)}
                    if len(value_set) < len(grouped[key]):
                        nan_keys.add(key) # Add key which had NaN
                    if key in data_dict:
                        data_dict[key].update(value_set)
                    else:
                        data_dict[key] = value_set
        with open(nan_file,"w") as nf:
            for key in nan_keys:
                nf.write(key + '\n')
        if nan_keys:
            print("NaN values were found while constructing dictionary you blast db or taxonomy db might be out of date." + "\n"
                  + f"Keys associated with NaN values were saved to: {nan_file} analysis can proceed without those values.")
        else:
            print("No NaN values detected")
        return data_dict
    @staticmethod
    def read_sequence_ids(directory_path):
        sequence_ids_dict = {}
        special_dict = {}

        directory = Path(directory_path)
        for file_path in directory.glob('unique_*'):
            with file_path.open('r') as file:
                sequence_ids = [line.split()[0] for line in file if line.strip()]
            
            if sequence_ids:
                if file_path.name == 'unique_taxids_cluster_sequences.txt' or file_path.name =="unique_taxids_filtered_sequences_singles.txt":
                    special_dict['special'] = {seq_id.split('_')[-1] for seq_id in sequence_ids}
                else:
                    sequence_ids_dict[file_path.stem] = {seq_id.split('_')[-1] for seq_id in sequence_ids}
            else:
                print(f'The file {file_path.name} is empty.')
        return sequence_ids_dict,special_dict
    @staticmethod
    def filter_fasta_file_dict(fasta_file_path,sequence_ids_dict,special_dict,output_directory):
        combined_ids = set()
        for ids in sequence_ids_dict.values():
            combined_ids.update(ids)
        special_ids = special_dict.get('special',set())

        fasta_file_path = Path(fasta_file_path)
        output_path = Path(output_directory) / (fasta_file_path.stem + "_filtered.fasta")
        with fasta_file_path.open('r') as fasta, output_path.open('w') as output:
            write_sequence = False
            for line in fasta:
                if line.startswith('>'):
                    seqID_last_number = line.split('>')[1].split('_')[-1].strip()
                    prepend_cluster = seqID_last_number in special_ids
                    if seqID_last_number in combined_ids or prepend_cluster:
                        output_line = f">CLUSTER_{line[1:]}" if prepend_cluster else line
                        write_sequence = True
                    else: write_sequence = False
                if write_sequence:
                    output.write(output_line if line.startswith('>') else line)
        print(f"Filtered FASTA file has been written to {output_path}")
        return output_path
    @staticmethod
    def create_fasta_file_with_full_TRS_sequences(df,filtered_df,results_directory,incrementent_indices):
        """
        Creates fasta file containing full sequences using filtered dataframe and the original one
        In addition assigns new unique ids to each sequence that contain information about flanking TRS
        """
        ncbi_ids = df["GENOME"].unique().tolist() # maybe move it into arguments ?
        organism_map = SequenceProcessor.fetch_organism_names(ncbi_ids)
        filtered_df['Taxonomic Name'] = df["GENOME"].map(organism_map)
        indices = incrementent_indices
        filtered_df['L_R_id'] = filtered_df['Taxonomic Name'] + '_L' + filtered_df['L-No'].astype(str) + '_R' + filtered_df['R-No'].astype(str)
        full_sequences_path = os.path.join(results_directory,"final_output")
        full_sequences_fasta = os.path.join(full_sequences_path,"full_sequences.fasta")
        with open(full_sequences_fasta,'w') as file:
            for index, row in filtered_df.iterrows():
                file.write(f'>{row["L_R_id"]}_{index}\n')
                file.write(f'{row[">SEQ"]}\n')
    @staticmethod
    def parse_fasta_indices(fasta_file):
        """
        Parse fasta file and extract the part after the last '_' from each sequence id. Increment those values by 1 and return them
        Those incremented values are the REAL positions of found sequences in csv file contataining all sequences.
        Remove the incrementing part when the 1 indexed ids are no longer used
        """
        incremented_indices = []
        indices = []
        with open(fasta_file,'r') as file:
            for line in file:
                if line.startswith('>'):
                    index = line.strip().split('_')[-1]
                    indices.append(index)
                    incremented_index = int(index) - 1
                    incremented_indices.append(incremented_index)
        incremented_indices = list(set(incremented_indices)) #Stupid,unnecessary but later code fixes it 
        return incremented_indices
    @staticmethod
    def extract_rows_from_csv(df,row_indices):
        """Filters the DataFrame using indices corresponding to rows in dataframe and sorts them using those indexes"""
        filtered_df = df.iloc[row_indices].sort_index()
        return filtered_df
    @staticmethod
    def extract_full_TRS_sequences(df,fasta_file,results_directory):
        """Executes the above functions"""
        #SequenceProcessor.validate_and_set_email() # Debug purposeses
        incremented_indices = BLASTProcessor.parse_fasta_indices(fasta_file)
        filtered_df = BLASTProcessor.extract_rows_from_csv(df=df,row_indices=incremented_indices)
        BLASTProcessor.create_fasta_file_with_full_TRS_sequences(df,filtered_df,results_directory,incrementent_indices=incremented_indices)
