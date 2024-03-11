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
import matplotlib.pyplot as plt
from multiprocessing import Pool
import re
import subprocess
import sys
import fnmatch
from tqdm import tqdm

class SeqProcessor:
    @staticmethod
    def ensure_directory_exists(directory_path):
        """
        Ensure that the specified directory exists. If it does not exist, it is created.
        
        Parameters:
        - directory_path (str): The path of the directory to check and potentially create.
        """
        try:
            if not os.path.exists(directory_path):
                os.makedirs(directory_path)
                print(f"Directory {directory_path} created successfully.")
            else: 
                print(f"Directory {directory_path} already exists.")
        except PermissionError:
            print(f"Permission denied: Unable to create {directory_path}.")
        except FileNotFoundError:
            print(f"One or more intermediate directories in the path do not exist")
        except Exception as e:
            print(f"An error occured while creating {directory_path} : {e}")
    
    @staticmethod
    def calculate_sequence_lengths(tmin):
        """
        Calculate maximum lengths for left and right sequences based on a minimum total length.
        
        Parameters:
        - tmin (int): The minimum total length for both sequences combined.
        
        Returns:
        Tuple of (l_chars_max, r_chars_max) representing the maximum lengths for left and right sequences.
        """
        try:
            l_chars_max = tmin // 2
            r_chars_max = tmin // 2
            print(f"The maximum available length for SEQ_L and SEQ_R is {l_chars_max} and {r_chars_max} respectively")
            return l_chars_max, r_chars_max
        except Exception as e:
            print(f"An error has occured: {e}")
    
    @staticmethod
    def adjust_input_to_range(prompt_message, max_val):
        """
        Promps the user for input and adjusts it to ensure it does not exceed the maximum allowed value.
        
        Parameters:
        - prompt_message (str): The message that user will see.
        - max_val (int): The maximum allowed value.
        
        Returns:
        The adjusted value, ensuring it does not exceed the maximum value.
        """
        while True:
            try:
                user_input = int(input(prompt_message))
                if user_input > max_val:
                    print(f"Your input was adjusted to {max_val}.")
                    return max_val
                return user_input
            except ValueError:
                print("Please ensure that you provide a valid integer")
            except Exception as e:
                print(f"An unexpected error has occured: {e}")
    @staticmethod
    def extract_sequences(combined_results, l_chars, r_chars):
        """
        Extracts left and right sequences from the combined results based on specified lengths.
        
        Parameters:
        - combined_results (DataFrame): The DataFrame containing the sequences.
        - l_chars (int): The length of the sequence to extract from the start.
        - r_chars (int): The length of the sequence to extract from the end.
        
        Returns:
        The modified DataFrame with new columns for the left and right sequences.
        """
        try:
            combined_results['SEQ_L'] = combined_results['>SEQ'].str.slice(0, l_chars)
            combined_results['SEQ_R'] = combined_results['>SEQ'].str.slice(-r_chars)
            print("Dataframe filtered successfully")
            return combined_results
        except Exception as e:
            print(f"An error has occured: {e}!")

    @staticmethod
    def validate_email(email):
        """
        Validate email format.

        Args:
            email(str): Email adress to validate.
        
        Returns:
            bool: True if the email format is valid, False otherwise
        """
        try:
            return bool(re.match(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$", email))
        except Exception as e:
            print(f"An error occured during email validation:{e}. Please use standard email formats.")
            return False       
    @staticmethod
    def set_user_email():
        """
        Prompt the user for their email adress and set it for Entrez access.
        """
        while True:
            try:
                print("To further protect your email adress please consider using a disposable or temporary one")
                user_email = input("Enter your email adress to access Entrez:")
                if SeqProcessor.validate_email(user_email):
                    Entrez.email = user_email
                    print("Email adress set successfully")
                    break
                else:
                    print("Invalid email format. Please enter a valid email adress.")
            except Exception as e:
                print(f"An unexpected error has occured: {e}. Please try again.")
    @staticmethod
    def fetch_organism_names(ncbi_ids):
        """
        Fetch organism names from NCBI based on NCBI IDs.

        Args:
            ncbi_ids (list): List of NCBI IDs.

        Returns:
            dict: Dictionary mapping NCBI IDs to organism names.
        """
        organism_map = {}
        for ncbi_id in ncbi_ids:
            try:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb",
                                       retmode ="text")
                record = handle.read()
                for line in record.splitlines():
                    if line.startswith("  ORGANISM"):
                        print(f"Fetching organism name for id: {ncbi_id}.")
                        organism = line.split("ORGANISM")[1].strip().replace(" ","_")
                        organism_map[ncbi_id] = organism
                        print(f"Retrieved organism name for id: {ncbi_id} successfuly.")     
            except Exception as e:
                print(f"Error fetching organism name for NCBI ID {ncbi_id}: {e}")
        print(f"All organism names successfully retrieved")
        return organism_map
    @staticmethod
    def rename_sequences(input_file, output_file):
         """
        Renames sequences in a FASTA file.

        Parameters:
            input_file (str): Path to the input FASTA file.
            output_file (str): Path to the output FASTA file with renamed sequences.
        """
         try:
            with open(input_file,"r") as input_handle, open(output_file,"w") as output_handle:
                # Initialize pair counter
                pair_index = 1
                for record in SeqIO.parse(input_handle,"fasta"):
                    orignal_id = record.id
                    new_id = f"{orignal_id}_{pair_index}"
                    record.id = new_id
                    record.description = ""
                    SeqIO.write(record,output_handle,"fasta")
                    if 'R' in orignal_id:
                        pair_index += 1
         except FileNotFoundError:
             print(f"Error: The file {input_file} was not found!")
         except PermissionError:
             print(f"Permission denied when accessing the files!")
         except IOError as e:
             print(f"An I/O error occured: {e}")
         except Exception as e:
             print(f"An unexpected error has occured: {e}")
         print(f"Seqeuences successfully renamed and saved at {output_file}!")
    @staticmethod
    def adjust_word_length(c):
        """
        Dynamically adjusts 'n' based on sequence identity threshold 'c'. 
        """
        try:
            threshold_n_values = [
                (0.95,10),
                (0.9,8),
                (0.88,7),
                (0.85,6),
                (0.8,5),
                (0.75,4),
            ]
            if c==1.0:
                return 11
            for lower_bound, n_value in threshold_n_values:
                if c >= lower_bound:
                    print(f"Wordsize was set to {n_value}")
                    return n_value
            return None
        except TypeError as e:
            print(f"Error:{e}")
        except ValueError as e:
            print(f"Error:{e}")
        except Exception as e:
            print(f"An unexpected error has occured: {e}")
    @staticmethod
    def run_cdhit(cdhit_path, input_file, output_file, c=None, d=0, m="t", g=0, G=1, sc=0, results_directory="results"):
        """
        Run cdhit directly from python adjusting the 'n' based on user provided c value.
        """
        if c is None:
            while True:
                try:
                    c = float(input("Enter the sequence identity threshold (between 0.75 and 1.0): "))
                    if 0.75 <= c <= 1.0:
                        break  # Valid input; break from the loop
                    else:
                        print("Error: The sequence identity threshold must be between 0.75 and 1.0.")
                except ValueError:
                    print("Invalid input. Please enter a numerical value.")
        
        n = SeqProcessor.adjust_word_length(c)  # Correctly call the class method

        if SeqProcessor.cdhit_path_global and os.path.exists(os.path.join(SeqProcessor.cdhit_path_global, "cd-hit")):
            cdhit_path = SeqProcessor.cdhit_path_global
        elif not cdhit_path or not os.path.exists(os.path.join(cdhit_path, "cd-hit")):
            cdhit_path = input("Enter the directory where CD-HIT is located: ")
            if not os.path.exists(os.path.join(cdhit_path, "cd-hit")):
                raise FileNotFoundError("CD-HIT executable not found in the specified directory.")
        SeqProcessor.cdhit_path_global = cdhit_path  # Correctly update the class variable
        
        SeqProcessor.ensure_directory_exists(results_directory)  # Correctly use the class method

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
    @staticmethod
    def run_cdhit_new(cdhit_path, input_file, output_file, c=None, d=0, m="t", g=0, G=1, sc=1, results_directory="results"):
        """
        """
        if c is None:
            while True:
                try:
                    c = float(input("Enter the sequence identity threshold between 0.75 and 1: "))
                    if 0.75 <= c <= 1.0:
                        break
                    else:
                        print("Error the value is out of range")
                except ValueError:
                    print("Invalid input. Please enter numerical value")
        
        n = SeqProcessor.adjust_word_length(c)

        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file {input_file} not found")
        
        if not os.path.isdir(results_directory):
            SeqProcessor.ensure_directory_exists(results_directory)
        
        if cdhit_path is None:
            cdhit_path = input("Enter the directory where CD-HIT is located: ")
        if not cdhit_path or not os.path.exists(os.path.join(cdhit_path,"cd-hit-est")):
            raise FileNotFoundError(f"CD-HIT executable was not found in specified directory:{cdhit_path}")
        
        cmd = [
            os.path.join(cdhit_path,"cd-hit-est"),
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

        try:
            subprocess.run(cmd,check=True)
            print(f"CD-HIT executed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"CD-HIT command failed: {e}")
        except OSError as e:
            print("Execution failed:{e}")
        
        try:
            new_results_directory = f"{results_directory}_c{c}"
            if not os.path.exists(new_results_directory):
                os.rename(results_directory, new_results_directory)
                print(f"The {results_directory} will be renamed to {new_results_directory}")
                results_directory = new_results_directory
                print(f"The results directory has been renamed to: {results_directory}")
            else:
                print(f"The directory {new_results_directory}, already exists!")
        except Exception as e:
            print(f"An error has occured while managing the results directory: {e}")

        return results_directory

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
                for line in lines:
                    sequence_id = line.split(">")[1].split("...")[0].strip()
                    outfile.write(">" + sequence_id + "\n")
        except FileNotFoundError:
            print(f"Error: file {input_file} does not exist!")
        except PermissionError:
            print(f"Error: Access to {input_file} denied")
        except Exception as e:
            print(f"An unexpected error has occured!")
        print(f"Overwriting file {input_file} with cleaned entries")
    @staticmethod
    def read_fasta_ids(filename):
        """
        Read and return a set of FASTA IDs from a given file.

        Parameters:
        - filename (str): Path to the FASTA file.

        Returns:
        - set: A set of FASTA IDs found in the file.
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
            print(f"Error: Permision denied when accessing {filename}")
        print("Fasta IDs processed successfully")
        return fasta_ids

    @staticmethod
    def filter_fasta_file(input_file, output_file, fasta_ids_to_remove, chunk_size=1024*1024):
        """
        Filter entries from a FASTA file, removing sequences with IDs in the given set.

        Parameters:
        - input_file (str): Path to the input FASTA file.
        - output_file (str): Path where the filtered FASTA file will be saved.
        - fasta_ids_to_remove (set): A set of FASTA IDs to be removed from the input file.
        - chunk_size (int): Size of the chunk to read at a time (in bytes). Default is 1MB.
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
        except PermissionError:
            print(f"Error: Permission denied when acessing the file {input_file} or {output_file}")
        except Exception as e:
            print(f"An error occurred: {e}")
    @staticmethod
    def filter_fasta_file_clusters(input_file, output_file, fasta_ids_to_include):
        """
        Filters a FASTA file to only include sequences with IDs present in the given set.

        Parameters:
        - input_file (str): Path to the input FASTA file.
        - output_file (str): Path where the filtered FASTA file will be saved.
        - fasta_ids_to_include (set): A set of FASTA IDs that should be included in the output file.
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
        except PermissionError:
            print(f"Error: Permission denied when accessing the file {input_file} or {output_file}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
    @staticmethod
    def find_directory_by_name(directory_name):
        """
        Searches for directories with a specified name in common locations on the user's computer 
        and allows the user to select one if multiple are found.

        Params: 
        -directory_name (str): The name of the directory to search for.

        Returns:
        -str : The selected path to directory

        Currently looks for all folders matching the name on the computer might need to think about 
        limiting it to a cwd or something 
        """
        common_locations = [
            os.path.expanduser('~'),
            os.path.join(os.path.expanduser('~'),'Documents'),
            os.path.join(os.path.expanduser('~'), 'Desktop')
        ]

        found_directories = []
        for location in common_locations:
            try:
                for root,dirs, _ in os.walk(location):
                    if directory_name in dirs:
                        found_directories.append(os.path.join(root,directory_name))
            except PermissionError as e:
                print(f"Permission denied when acessing {location}. Error: {e}")
            except Exception as e:
                print(f"An unexprected error has occured: {e}")
        
        #If no directories were found
        if not found_directories:
            print("No directories found check your input")
            return None

        #If only one directory was found return it
        if len(found_directories) == 1:
            return found_directories[0]
        
        #If multiple directories with similar names were found prompt user for input
        print(f"Mulitple {directory_name} directories found. Please select one:")
        for i,directory in enumerate(found_directories,start=1):
            print(f"{i}. {directory}")
        
        while True:
            selection = input("Enter the number of the directory you want to use (or exit by typing 'exit'): ")
            if selection.lower() == 'exit':
                print("Exiting.")
                return None
            
            try:
                selected_index = int(selection) - 1 #Python is zero-based when it comes to indexes
                if 0 <= selected_index < len(found_directories):
                    return found_directories[selected_index]
                else:
                    print("Invalid selection. Try again.")
            except ValueError:
                print("Invalid input. Please enter a number")
    @staticmethod
    def search_files(search_path,file_name):
        found = []
        for root, _, files in os.walk(search_path):
            for file in files:
                if fnmatch.fnmatch(file,file_name):
                    found.append(os.path.join(root,file))
        return found
    @staticmethod
    def prompt_user_selection(found_files,file_name,auto=False):
        if not found_files:
            print("No files/directories with specified names were found.")
            return None
        
        if auto or len(found_files) == 1:
            return found_files[0]
        
        print(f"Multiple {file_name} files/directories found. Please select one:")
        for i, file in enumerate(found_files,start=1):
            print(f"{i}. {file}")
        
        while True:
            selection = input("Enter the number of the file/directory you want to use or type 'exit' to exit:")
            if selection.lower() == "exit":
                print("Exiting.")
                return None

            try:
                selected_index = int(selection) - 1
                if 0 <= selected_index < len(found_files):
                    return found_files[selected_index]
                else:
                    print(f"Invalid selection. Try again.")
            except ValueError:
                print("Invalid input. Please enter a number")
    @staticmethod
    def find_file_by_name(file_name, auto=False, folder=None):
        """
        Searches for files with a specified name. If neither 'auto' nor 'folder' is specified,
        searches common locations and prompts the user to select a file if multiple are found.
        
        Args:
            file_name (str): The name of the file to search for.
            auto (bool): If True, automatically returns the path of the first file found.
            folder (str): Optional. If provided, limits the search to this specified path.

        Returns:
            str or None: Path of the selected or found file, or None if not found or exited.
        """
        if folder:
            #If folder is specified limit the search to directory
            found_files = SeqProcessor.search_files(folder,file_name)
            return SeqProcessor.prompt_user_selection(found_files,file_name,auto=auto)
        elif auto:
            #If auto is set to True but no folder is specified, search the common locations for first match
            common_locations = [
                os.path.expanduser('~'),
                os.path.join(os.path.expanduser('~'), 'Documents'),
                os.path.join(os.path.expanduser('~'), 'Desktop')
            ]
            for location in common_locations:
                found_files = SeqProcessor.search_files(location,file_name)
                if found_files:
                    return found_files[0]
                else:
                    print("No files with specified name found in common locations")
        else:
            common_locations = [
            os.path.expanduser('~'),
            os.path.join(os.path.expanduser('~'), 'Documents'),
            os.path.join(os.path.expanduser('~'), 'Desktop')
            ]
            all_found = []
            for location in common_locations:
                all_found.extend(SeqProcessor.search_files(location,file_name))
            return SeqProcessor.prompt_user_selection(all_found,file_name,auto=auto)    
            
    @staticmethod
    def filter_and_overwrite_blast_file(file_path):
        """
        Processes BLAST output file to remove duplicate pairs of sequence IDs and accession numbers,
        keeps the first occurence of each unique pair along with the entire line of data.

        This function reads a file specified by `file_path`, identifies the first occurrence of unique 
        sequence ID and accession number pairs, and overwrites the original file with those lines.
    
        Parameters:
        - file_path (str): The path to the BLAST output file to be processed.
        """
        unique_pairs = set()
        lines_to_keep = []

        try:
            with open(file_path,"r") as file:
                for line in file:
                    try:
                        columns = line.strip().split('\t')
                        if len(columns) >= 2:
                            pair = (columns[0],columns[1])
                            if pair not in unique_pairs:
                                unique_pairs.add(pair)
                                lines_to_keep.append(line)
                    except Exception as e:
                        print(f"An unexpected error has occured while processing a line in {file_path}: {e}")
                print(f"Searching for unique sequence ID - acession pairs in {file_path}....")
        except PermissionError as e:
            print(f"Permission denied while trying to access {file_path}: {e}")
            return
        except IOError as e:
            print(f"An error has occured while opening or reading {file_path}: {e}")
            return
        
        try:
            with open(file_path,"w") as file:
                file.writelines(lines_to_keep)
        except PermissionError as e:
            print(f"Permission denied while to access {file_path} for writing: {e}")
        except IOError as e:
            print(f"An error occured while writing to {file_path}: {e}")
    @staticmethod
    def convert_to_txt(directory_path):
        """
        Identifies non-txt files in directory_path and converts them to tab separated files with .txt extension,
        after user confirmation if found files are valid

        Parameters:
        - directory_path(str): The path to the directory where files will be verified and potentially renamed.
        """
        try:
            non_txt_files = [f for f in os.listdir(directory_path) if not f.endswith(".txt")]
        except FileNotFoundError:
            print("The specified directory does not exist.")
            return
        except PermissionError:
            print(f"Permission denied to {directory_path}.")
            return
        except Exception as e:
            print(f"An unexpected error has occured: {e}")

        if non_txt_files:
            print(f"No .txt files found in {directory_path}")
            for file in non_txt_files:
                print(f"Found file: {file} lacks extension.")
            
            print("Directory should contain only blast results!")
            user_input = input(f"Do you want to rename these files to include .txt extension? (y/n): ")
            if user_input.lower() == "y":
                for file in non_txt_files:
                        orignal_path = os.path.join(directory_path,file)
                        new_path = f"{orignal_path}.txt"
                        try:
                            os.rename(orignal_path,new_path)
                            print(f"Renamed {file} to {file}.txt")
                        except FileNotFoundError:
                            print(f"Failed to rename {file}: File not found.")
                        except PermissionError:
                            print(f"Failed to rename {file} : Permission denied.")
                        except Exception as e:
                            print("Failed to rename {file}: {e}")
            elif user_input.lower() == "n":
                print("No files were renamed.")
            else:
                print("Invalid input. Please type 'y' for yes or 'n' for no.")
        else:
            print("No non-txt files were found.")
    @staticmethod
    def process_blast_files_in_directory(directory_path):
        """
        Processes all BLAST output files in a specified directory to remove duplicate pairs,
        keeping the first occurrence of each unique pair along with the entire line of data.
        
        This function searches for all files in `directory_path` ending with '.txt', and for each file,
        it calls `filter_and_overwrite_blast_file` to process and overwrite it with a filtered version
        containing only unique sequence ID and accession number pairs.
        
        Parameters:
        - directory_path (str): The path to the directory containing BLAST output files.
        """
        SeqProcessor.convert_to_txt(directory_path)

        txt_files = [f for f in os.listdir(directory_path) if f.endswith(".txt")]
    
        if not txt_files:
            print("No .txt files found in the directory after attempting to rename. Exiting program.")
            sys.exit("Exiting program because no .txt files were found.")

        for filename in txt_files:
            try:
                file_path = os.path.join(directory_path,filename)
                SeqProcessor.filter_and_overwrite_blast_file(file_path)
                print(f"Processed {filename}")
            except Exception as e:
                print(f"An unexpected error has occured: {e}!")
    @staticmethod
    def collect_accessions_from_blast_files(folder_path):
        """Collects unique acessions from BLAST files in the specified folder"""
        unique_accessions = set()
        try:
            for filename in os.listdir(folder_path):
                if filename.endswith(".txt"):
                    file_unique_accessions = set()
                    file_path = os.path.join(folder_path,filename)
                try:
                    with open(file_path,"r") as file:
                        for line in file:
                            try:
                                columns = line.strip().split('\t')
                                if len(columns) >= 2:
                                    accession = columns[1]
                                    file_unique_accessions.add(accession)
                                    unique_accessions.add(accession)
                            except ValueError as e:
                                print(f"Error processing line in {filename}: {e}")
                    print(f"{filename}: Found {len(file_unique_accessions)} unique accessions.")
                except IOError as e:
                    print(f"Error opening or reading file {filename}: {e}")
        except Exception as e:
            print(f"An unexpected error has occured")
        
        print(f"Total unique accessions across all files: {len(unique_accessions)}")
        return unique_accessions
    @staticmethod
    def filter_taxonomy_file(taxonomy_file,accessions,chunksize=50000):
        """
        Filters the taxonomy file in chunks, keeping rows matching the accessions from BLAST files. 
        Informs of overall chunk progress
        Params:
        - taxonomy_file (str) - path to processed taxonomy database
        - accessions (lst) - list of acessions in our files
        - chunksize (int) - size of each chunk  
        """
        try:
            filtered_rows = []
            cols_to_use = [1,2]
            chunk_count = 0
            estimated_chunks = 6548

            for chunk in tqdm(pd.read_csv(taxonomy_file, sep='\t', usecols=cols_to_use, header=None, chunksize=chunksize), total=estimated_chunks):
                #chunk_count += 1
                filtered_chunk = chunk[chunk.iloc[:, 0].isin(accessions)]
                progress_percentage = (chunk_count / estimated_chunks) * 100
                #print(f"Progress: {progress_percentage:.2f}% processed")
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
            print(f"An unexpected error has occured: {e}")
    @staticmethod
    def map_accession_to_taxid(accession): #deprecated
        """Map accession to TaxID using the optimized lookup dictionary"""
        return accession.get(accession,'') # Empty string if accession not found
    @staticmethod
    def fetch_organism_taxids(ncbi_ids):
        """
        Fetch organism TaxIDs from NCBI based on NCBI IDs, using organism names as a fallback.
        """
        organism_taxid_map = {}
        organism_names = SeqProcessor.fetch_organism_names(ncbi_ids)

        for ncbi_id in ncbi_ids:
            try:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
                record = handle.read()
                handle.close()
                taxid = None
                for line in record.splitlines():
                    if line.strip().startswith("/db_xref=\"taxon:"):
                        taxid = line.split(":")[1].replace("\"", "").strip()
                        break
                if taxid:
                    organism_taxid_map[ncbi_id] = taxid
                    print(f"Retrieved TaxID {taxid} for NCBI ID {ncbi_id} successfully.")
                else:
                    print(f"TaxID not found for {ncbi_id}")
            except Exception as e:
                print(f"Error fetching TaxID for NCBI ID {ncbi_id}: {e}")

        # Convert keys in organism_taxid_map from NCBI IDs to organism names
        organism_taxid_map_with_names = {}
        for ncbi_id, taxid in organism_taxid_map.items():
            organism_name = organism_names.get(ncbi_id, "Unknown_Organism")
            organism_taxid_map_with_names[organism_name] = taxid

        print("All organism TaxIDs successfully retrieved and mapped to organism names.")
        return organism_taxid_map_with_names
    @staticmethod
    def filter_key_parts(key):
        """
        Filters parts of the key based on the following rules:
        - Drops the part if it contains a number
        - Drops the part if it starts with a capital letter and is not in the first part 
        """
        parts = key.split("_")
        filtered_parts = [parts[0]] #Always make sure that first part is already included

        for part in parts[1:]:
            if not re.search(r'^[A-Z]',part) and not re.search(r'\d',part):
                filtered_parts.append(part)
        
        return '_'.join(filtered_parts)
    @staticmethod
    def interact_and_update_dict(dictionary):
        """
        Interacts with user to possibly add more TaxIDs and organism names to species_info dictionary

        Args:
            dictionary(dict): The exisitng dictionary mapping organism names to SETS of TaxIDs.
        """
        user_decision = input("Would you like to provide additonal NCBI ID to include in dictionary? (y/n)")

        if user_decision.lower() == 'y':
            user_input = input("Please enter the NCBI IDs as comma separated list or a single value: ")
            ids = [id.strip() for id in user_input.split(",")]
            print("Stripping the list completed.")
            new_organism_ids = SeqProcessor.fetch_organism_taxids(ids)
            for organism_name, taxid_set in new_organism_ids.items():
                if organism_name in dictionary:
                    dictionary[organism_name].update(taxid_set)
                else:
                    dictionary[organism_name] = taxid_set
        
        elif user_decision.lower() == 'n':
            print("No additional NCBI Ids will be added.")
        else:
            print("Invalid response. Nothing will be added...")
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
    def remove_keys_with_multiple_values(dict):
        """
        Self explanatory 

        Args
        -dict (dict) - dictionary to remove multiple values from
        """
        unique_dict = {}

        #Iterate over items in the dict
        for key, value in dict.items():
            #Check if the value associted with key is a list and if it has more than one entry:
            if isinstance(value,list) and len(value) == 1:
                #keep only the first one
                unique_dict[key] = value[0]
            elif not isinstance(value,list):
                unique_dict[key] = value
        return unique_dict
    @staticmethod
    def process_files_with_filter(directory,filtered_dict):
        output_directory = os.path.join(directory,"unique_sequences")
        SeqProcessor.ensure_directory_exists(output_directory)

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
    def process_taxid_input(input_str):
        """
        Process the input string which could be single taxid, comma separated list of taxids or taxid file with one taxid per row
        """
        taxids = set()
        if input_str.endswith('.txt'):
            input_str = SeqProcessor.find_file_by_name(input_str)
            try:
                with open(input_str,'r') as file:
                    line_number = 0 
                    for line in file:
                        line_number +=1
                        taxid = line.strip()
                        if taxid.isdigit():
                            taxids.add(int(taxid))
                        else:
                            print(f"Warning: Line {line_number} in file {input_str} is badly formatted: {taxid}")
            except FileExistsError:
                print(f"File {input_str} not found")
        else:
            #Process a single taxid or comma separated list of taxids
            for taxid in input_str.split(','):
                taxid = taxid.strip()
                if taxid.isdigit():
                    taxids.add((taxid))
                else:
                    print(f"Warning: Badly formatted input detected ('{taxid}').")
        return taxids
    @staticmethod
    def ask_for_exceptions():
        exceptions = set()
        while True:
            user_input = input("Would you like to add exceptions to further filtering steps(ex. top-level species) (y/n?))")
            if user_input.lower() == 'y':
                print("Enter a taxid, a list of comma-separated taxids, or a file path with one taxid per row. Note: If entering a file path, ensure it ends with '.txt'.")
                taxid_input = input("Enter your choice here: ")
                exceptions = SeqProcessor.process_taxid_input(taxid_input)
                break
            elif user_input.lower() == 'n':
                print("No exceptions will be added.")
                break
            else:
                print("Invalid input. Please enter 'y' for yes and 'n' for no")
        return exceptions
    @staticmethod
    def filter_with_exceptions(results_dict,species_info,exceptions):
        '''
        Filters the results dictionary based on the presence of taxids in the species infor dictionary,
        taking into account exceptions.

        Args:
            results_dict(dict) - dictionary with keys as sequence ids and sets of values as taxids
            species_info(dict) - dictionary with species names as keys and sets of taxids as values
            exceptions(set) - set of taxids that are exceptions
        
        Returns:
            dict: Filtered results dictionary
        '''
        filtered_keys = {}
        all_species_taxids = set.union(*species_info.values()) # combine all taxids from species info for easier lookup

        for key, values in results_dict.items():
            non_exception_values = values - exceptions # Remove exceptions 
            if non_exception_values: #If there are none left
                #Keep keys that have exactly one value that matches the species info
                if len(non_exception_values) == 1 and non_exception_values & all_species_taxids:
                    filtered_keys[key] = non_exception_values
            else: #If all values are exceptions 
                #Check if ANY of the exceptions are in the species_info dictionary(indicating we should keep them)
                if values & all_species_taxids:
                    filtered_keys[key] = values
        return filtered_keys
    @staticmethod
    def count_fasta_sequences(file_path):
        try:
            with open(file_path,"r") as fasta_file:
                return sum(1 for line in fasta_file if line.startswith('>'))
        except FileNotFoundError:
            print(f"File {file_path} not found.")
            return None
    @staticmethod
    def count_clusters(file_path):
        '''
        Counts the total number of clusters and the number of clusters with more than 2 entires in
        CD-HIT cluster file
        '''
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

                    #Check for the last cluster
                    if current_cluster_size > 2:
                        clusters_with_more_than_two_entries +=1
        
        except FileNotFoundError:
            print(f"File {file_path} not found")
            return None,None
        
        return total_clusters,clusters_with_more_than_two_entries
    
    @staticmethod
    def count_lines_in_file(file_path):
        try:
            with open(file_path,"r") as file:
                return sum(1 for _ in file)
        except FileNotFoundError:
            print(f"File {file_path} not found.")
            return None
    
    @staticmethod
    def write_summary(summary_file_path, results_directory):
        """
        Writes a summary including fasta files statistics variables defined etc.
        """
        optional_variables = {}
        for var in ['tmin','tmax','mode','exceptions']:
            try:
                optional_variables[var] = eval(var)
            except NameError:
                pass
        
        expected_files = [
            'combined_sequences.fasta',
            'cluster_sequences_combined_unique.fasta',
            'filtered_sequences_combined_unique.fasta',
            'combined_sequences_unique_cdhit.clstr',
            'filtered_sequences_combined_unique_blastn_out.txt',
            'cluster_sequences_combined_unique_blastn_out.txt',
            'taxids_cluster_sequences_combined_unique_blastn_out.txt',
            'taxids_filtered_sequences_combined_unique_blastn_out.txt'
        ]

        found_files_info = []

        for file_name in expected_files:
            file_path = SeqProcessor.find_file_by_name(file_name,results_directory)
            if file_path:
                if file_name.endswith('.fasta'):
                    num_sequences = SeqProcessor.count_fasta_sequences(file_path)
                    found_files_info.append(f"{file_name} contains {num_sequences} sequences.")
                elif file_name.endswith('.clst'):
                    total_clusters, clusters_with_more = SeqProcessor.count_clusters(file_path)
                    found_files_info.append(f"CD-HIT Cluster file analysis for {file_name}: Total clusters {total_clusters}, Clusters with more than two entries {clusters_with_more}.")
                elif file_name.endswith('.txt'):
                    num_lines = SeqProcessor.count_lines_in_file(file_path)
                    found_files_info.append(f"{file_name} contains {num_lines} lines of data.")

        with open(summary_file_path, "w") as summary_file:
            # Write variables if present
            if optional_variables:
                summary_file.write("Execution Parameters:\n")
                for key in optional_variables.items():
                    summary_file.write(f"{key} : {value}\n")
                summary_file.write("\n")
            
            if found_files_info:
                summary_file.write("Found Files Information:\n")
                for info in found_files_info:
                    summary_file.write(f"{info}\n")
    @staticmethod
    def filter_and_overwrite_files_in_directory(directory_path, file_pattern="*.txt"):
        """
        Filters each file in a specified directory based on the last number in the identifiers of its lines,
        then overwrites each file with its filtered content. Targets files matching a given pattern.

        Args:
            directory_path (str): Path to the directory containing the files to process.
            file_pattern (str): Pattern of the file names to process (default is "*.txt").
        """

        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                with open(file_path, 'r') as file:
                    lines = file.readlines()

                number_count = {}
                filtered_entries = []

                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        number_count[number] = number_count.get(number, 0) + 1

                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit() and number_count[number] > 1:
                        filtered_entries.append(line.strip())

                with open(file_path, 'w') as file:
                    for entry in filtered_entries:
                        file.write(entry + '\n')
    @staticmethod
    def search_directories(search_path,directory_name):
        found_directories = []
        search_path = Path(search_path)
        for dir in search_path.glob('**/*'):
            if dir.is_dir() and fnmatch.fnmatch(dir.name,directory_name):
                found_directories.append(str(dir))
        return found_directories
    @staticmethod
    def find_directory_by_name_new(directory_name, auto= False, folder = None):
        '''
        This uses path to instead of os.walk like in previous one for now i'll leave it as an alternative function 
        Currently used only in part associated with blast
        Searches for directories with specified name.

        Args:
            directory_name(str) : The name of directory to search for
            auto(bool) : If true automatically select the first match
            folder(str): Provides the limits to the search (dir name to search in)
        '''
        if folder:
            folder = Path(folder)
            found_directories = SeqProcessor.search_directories(folder,directory_name)
        else:
            common_locations = [
                Path.home(),
                Path.home() / 'Desktop',
                Path.home() / 'Documents'
            ]
            all_found = []
            for location in common_locations:
                all_found.extend(SeqProcessor.search_directories(location,directory_name))
            found_directories = all_found

        return SeqProcessor.prompt_user_selection(found_directories,directory_name,auto=auto)   
    @staticmethod
    def read_sequence_ids(directory_path):
            sequence_ids_dict = {}
            special_dict = {}

            directory = Path(directory_path)
            for file_path in directory.glob('*blastn_out.txt'):
                with file_path.open('r') as file:
                    sequence_ids = [line.split()[0] for line in file if line.strip()]
                
                if sequence_ids:
                    if 'unique_taxids_cluster_sequences_combined_unique_blastn_out.txt' in file_path.name:
                        special_dict['special'] = {seq_id.split('_')[-1] for seq_id in sequence_ids}
                    else:
                        sequence_ids_dict[file_path.stem] = {seq_id.split('_')[-1] for seq_id in sequence_ids}
                else:
                    print(f'The file {file_path.name} is empty.')
            return sequence_ids_dict,special_dict
    @staticmethod
    def filter_fasta_file(fasta_file_path,sequence_ids_dict,special_dict,output_directory):
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
                        output_line = f">CLUSTER_{line[1:]}" if prepend_cluster else lin
                        write_sequence = True
                    else: write_sequence = False
                if write_sequence:
                    output.write(output_line if line.startswith('>') else line)
        print(f"Filtered FASTA file has been written to {output_path}")
