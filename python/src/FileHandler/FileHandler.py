import os 
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import fnmatch
import sys
from src.SequenceProcessor import SequenceProcessor

class FileHandler:
    
    @staticmethod
    def ensure_directory_exists(directory_path):
        '''
        Ensures that the specified directory exists if it does not creates it. 
        '''
        try:
            if not os.path.exists(directory_path):
                os.makedirs(directory_path)
            else:
                print(f"Directory {directory_path} already exists")
        except PermissionError:
            print(f"Premission denied: Unable to create {directory_path}")
        except FileNotFoundError:
            print(f"One or more intermediate directories in the path do not exists")
        except Exception as e:
            print(f"Unidentified error has occured while creating {directory_path}: {e}")
    
    @staticmethod
    def convert_to_txt(directory_path):
        '''
        Identifies non-txt files in the directory and converts them to tab-separated .txt files
        '''
        try:
            non_txt_files = [f for f in os.listdir(directory_path) if not f.endswith(".txt")]
        except FileNotFoundError:
            print("The specified directory does not exist.")
            return
        except PermissionError:
            print(f"Permission denied to {directory_path}.")
            return
        except Exception as e:
            print(f"An unexpected error has occurred: {e}")
        
        if non_txt_files:
            print(f"No .txt files found in {directory_path}")
            for file in non_txt_files:
                print(f"Found file: {file} lacks extension")

            print("Directory should contain only blast results!")
            for file in non_txt_files:
                original_path = os.path.join(directory_path,file)
                new_path = f"{original_path}.txt"
                try:
                    os.rename(original_path,new_path)
                    print("Renamed {file} to {file}.txt")
                except FileNotFoundError:
                    print(f"Failed to rename {file}: File not found.")
                except PermissionError:
                    print(f"Failed to rename {file}: Permission denied.")
                except Exception as e:
                    print(f"Failed to rename {file}: {e}")
        else:
            print("No non-txt files were found")

    @staticmethod
    def filter_and_overwrite_files_in_directory(directory_path, file_pattern= "*.txt"):
        """
        Filters each file in a specified directory based on the last number in the identifiers of its lines,
        then overwrites each file with its filtered content for entries with paired identifiers.
        Writes entries with unique identifiers to a separate file named '<original_filename>_singles.txt'.
        Targets files matching a given file pattern (default is "*.txt").
        """       
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename,file_pattern):
                file_path = os.path.join(directory_path,filename)
                single_entries_path = os.path.join(directory_path,f"{os.path.splitext(filename)[0]}_singles.txt")

                with open(file_path,'r') as file:
                    lines = file.readlines()

                number_count = {}
                filtered_entries = []
                non_unique_entries = []

                #First pass to count the occurences of each identifier
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        number_count[number] = number_count.get(number,0) + 1
                

                #Second pass to separate unique and non-unique identifiers
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        if number_count[number] > 1:
                            filtered_entries.append(line.strip())
                        else:
                            non_unique_entries.append(line.strip())
                
                # Overwrite the original file with paired identifiers
                with open(file_path, 'w') as file:
                    for entry in filtered_entries:
                        file.write(entry + '\n')

                # Write out single entries to a separate file
                with open(single_entries_path, 'w') as file:
                    for entry in non_unique_entries:
                        file.write(entry + '\n')
    
    @staticmethod
    def find_directory_by_name(directory_name, search_paths = None,auto=False):
        """
        Searches for directories with a specified name in provided locations and allows the user to select one if multiple are found.

        Params: 
        - directory_name (str): The name of the directory to search for.
        - search_paths (list of str): Optional. List of paths to search in. Defaults to common locations.
        - auto (bool): If True, automatically select the first found directory.
        
        Returns:
        - str : The selected path to the directory.
        """

        if search_paths is None:
            user_dir = os.path.expanduser('~')
            search_paths = [
                user_dir,
                os.path.join(user_dir,'Documents'),
                os.path.join(user_dir,'Desktop')
            ]
        
        found_directories = []

        for location in search_paths:
            try:
                for root,dirs, _ in os.walk(location):
                    if directory_name in dirs:
                        found_directories.append(os.path.join(root,directory_name))
            except PermissionError as e:
                print(f"Permission denied when accessing {location}. Error: {e}")
            except Exception as e:
                print(f"An unexpected error has occurred: {e}")

        if not found_directories:
            print("No directories found check your input.")
            return None
        
        if len(found_directories) == 1 or auto:
            return found_directories[0]
        
        print(f"Multiple {directory_name} directories found. Please select one:")
        for i, directory in enumerate(found_directories,start=1):
            print(f"{i}. {directory}")

        while True:
            selection = input("Enter the number of the directory you want to use(or type 'exit' to exit): ")
            if selection.lower() == 'exit':
                print('Exiting')
                return None
            
            try:
                selected_index = int(selection) - 1
                if 0 <= selected_index < len(found_directories):
                    return found_directories[selected_index]
                else:
                    print("Invalid selection. Try again")
            except ValueError:
                print("Invalid input. Please enter a number.")
    
    @staticmethod
    def read_fasta_ids(filename):
        """
        Read and return a set of FASTA IDs for a given file.
        """
        fasta_ids = set()

        def process_line(line):
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
        except PermissionError:
            print(f"Error: Permission denied when accessing {filename}")
        print("Fasta IDs processed successfully")
        return fasta_ids
    @staticmethod
    def parse_fasta_indices(fasta_file):
        """
        Parse fasta file and extract the part after the last '_' from each sequence id. 
        Increment those values by 1 and return them.
        """
        incremented_indices = []
        indices = []
        with open(fasta_file, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    index = line.strip().split('_')[-1]
                    indices.append(index)
                    incremented_index = int(index) - 1
                    incremented_indices.append(incremented_index)
        incremented_indices = list(set(incremented_indices))
        return incremented_indices
    @staticmethod
    def extract_rows_from_csv(df,row_indices):
        """
        Filters the DataFrame using indices corresponding to rows in dataframe and sorts them using those indexes.
        """
        filtered_df = df.iloc[row_indices].sort_index()
        return filtered_df
    @staticmethod
    def create_fasta_file_with_full_TRS_sequences(df,filtered_df, results_directory,incrementent_indices):
        """
        Creates fasta file containing full sequences using filtered dataframe and the original one.
        """
        ncbi_ids = df["GENOME"].unique().tolist()
        organism_map = SequenceProcessor.fetch_organism_names(ncbi_ids)
        filtered_df['Taxonomic Name'] = df["GENOME"].map(organism_map)
        indices = incrementent_indices
        filtered_df["L_R_id"] = filtered_df['Taxonomic Name'] + '_L' + filtered_df['L-No'].astype(str) + '_R' + filtered_df['R-No'].astype(str)
        full_sequences_path = os.path.join(results_directory, "final_output")
        full_sequences_fasta = os.path.join(full_sequences_path,"full_sequences.fasta")
        with open(full_sequences_fasta, 'w') as file:
            for index, row in filtered_df.iterrows():
                file.write(f'>{row["L_R_id"]}_{index}\n')
                file.write(f'{row[">SEQ"]}\n')
    @staticmethod
    def filter_and_overwrite_blast_file(directory_path):
        """
        Processes BLAST output files in the given directory to remove duplicate pairs of sequence IDs and accession numbers,
        keeping the first occurrence of each unique pair along with the entire line of data.
        """
        for filename in os.listdir(directory_path):
            if filename.endswith(".txt"):  # Process only .txt files
                file_path = os.path.join(directory_path, filename)
                unique_pairs = set()
                lines_to_keep = []

                try:
                    with open(file_path, 'r') as file:
                        for line in file:
                            try:
                                columns = line.strip().split('\t')
                                if len(columns) >= 2:
                                    pair = (columns[0], columns[1])
                                    if pair not in unique_pairs:
                                        unique_pairs.add(pair)
                                        lines_to_keep.append(line)
                            except Exception as e:
                                print(f"An unexpected error occurred while processing the line in {file_path}: {e}")
                    
                    # Write the filtered lines back to the file or a new file
                    with open(file_path, 'w') as file:
                        file.writelines(lines_to_keep)

                    print(f"Filtered unique sequence ID - accession pairs in {file_path}.")
                
                except PermissionError as e:
                    print(f"Permission denied while trying to access {file_path}: {e}")
                except IOError as e:
                    print(f"An error occurred while opening or reading {file_path}: {e}")
    @staticmethod
    def filter_fasta_file(input_file, output_file, fasta_ids_to_remove, chunk_size=1024*1024):
        FileHandler.filter_fasta(input_file, output_file, fasta_ids_to_remove, include=False, chunk_size=chunk_size)
    @staticmethod
    def filter_fasta_file_clusters(input_file, output_file, fasta_ids_to_include):
        FileHandler.filter_fasta(input_file, output_file, fasta_ids_to_include, include=True)
    @staticmethod
    def filter_fasta(input_file, output_file, ids_to_check, include=True, chunk_size=1024*1024):
        """
        Filter entries from a FASTA file based on IDs.

        Parameters:
        - input_file (str): Path to the input FASTA file.
        - output_file (str): Path where the filtered FASTA file will be saved.
        - ids_to_check (set): A set of FASTA IDs to be checked.
        - include (bool): If True, include sequences with IDs in ids_to_check. If False, exclude them.
        - chunk_size (int): Size of the chunk to read at a time (in bytes). Default is 1MB.
        """
        try:
            with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
                fasta_entry = []
                write_entry = include

                for line in f_in:
                    if line.startswith('>'):
                        if fasta_entry and write_entry:
                            f_out.writelines(fasta_entry)
                        fasta_entry = [line]
                        current_fasta_id = line.strip()[1:]
                        write_entry = (current_fasta_id in ids_to_check) if include else (current_fasta_id not in ids_to_check)
                    else:
                        fasta_entry.append(line)

                if fasta_entry and write_entry:
                    f_out.writelines(fasta_entry)
    
        except FileNotFoundError:
            print(f"Error: File {input_file} not found.")
        except PermissionError:
            print(f"Error: Permission denied when accessing the file {input_file} or {output_file}")
        except Exception as e:
            print(f"An error occurred: {e}")
    @staticmethod
    def search_files(search_path,file_name):
        found = []
        for root, _, files in os.walk(search_path):
            for file in files:
                if fnmatch.fnmatch(file,file_name):
                    found.append(os.path.join(root,file))
    @staticmethod
    def find_file_by_name(file_name, folder):
        """
        Searches for files with a specified name pattern within a given folder.

        Args:
            file_name (str): The name of the file to search for, supports wildcards.
            folder (str): The directory to search within.

        Returns:
            list: A list of file paths that match the search pattern, or an empty list if no files are found.
        """
        if not file_name or not isinstance(file_name, str):
            print("Invalid or empty file name provided")
            return []
        
        if not folder or not os.path.isdir(folder):
            print("Invalid or empty folder provided")
            return []
        
        found_files = []
        for root, dirs, files in os.walk(folder):
            for name in fnmatch.filter(files, file_name):
                found_files.append(os.path.join(root, name))
        
        return found_files
    @staticmethod
    def process_files_with_filter(directory,filtered_dict):
        output_directory = os.path.join(directory,"unique_sequences")
        FileHandler.ensure_directory_exists(output_directory)

        for filename in os.listdir(directory):
            if filename.endswith('.txt'):
                input_filepath = os.path.join(directory,filename)
                output_filepath = os.path.join(output_directory,f"unique_{filename}")

                df = pd.read_csv(input_filepath,sep='\t',header=None)

                #Create a boolean mask for rows where the first column matches a key in filtered_dict
                mask = df[0].isin(filtered_dict.keys())

                #Instead of creating a filtered dataframe we will directly alter the last column
                df.loc[mask,df.columns[-1]] = df.loc[mask,0].map(filtered_dict)

                #Now, filter the DataFrame to keep only the rows that match the mask
                filtered_df = df[mask]
                unique_df = filtered_df.drop_duplicates(subset=[0])
                unique_df.to_csv(output_filepath, sep='\t', index= False, header=False)
                print(f"Final output saved at {output_filepath}")
    # def process_files_with_filter(directory, filtered_dict, chunk_size=10000):
    #     output_directory = os.path.join(directory, "unique_sequences")
    #     FileHandler.ensure_directory_exists(output_directory)

    #     for filename in os.listdir(directory):
    #         if filename.endswith('.txt'):
    #             input_filepath = os.path.join(directory, filename)
    #             output_filepath = os.path.join(output_directory, f"unique_{filename}")

    #             # Set to accumulate unique rows based on the first column
    #             unique_rows = set()

    #             # Read the file in chunks
    #             reader = pd.read_csv(input_filepath, sep='\t', header=None, chunksize=chunk_size)
    #             for chunk in reader:
    #                 # Create a boolean mask for rows where the first column matches a key in filtered_dict
    #                 mask = chunk[0].isin(filtered_dict.keys())

    #                 # Directly alter the last column for the rows matching the mask
    #                 chunk.loc[mask, chunk.columns[-1]] = chunk.loc[mask, 0].map(filtered_dict)

    #                 # Filter the DataFrame to keep only the rows that match the mask
    #                 filtered_chunk = chunk[mask]

    #                 # Convert the DataFrame rows to tuples and add to the set of unique rows
    #                 unique_rows.update(filtered_chunk.itertuples(index=False, name=None))

    #             # Convert the set of unique rows back to a DataFrame
    #             unique_df = pd.DataFrame(list(unique_rows))

    #             # Write the unique DataFrame to the output file
    #             unique_df.to_csv(output_filepath, sep='\t', index=False, header=False)

    #             print(f"Final output saved at {output_filepath}")
    
