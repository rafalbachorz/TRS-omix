import os
import sys




class Stats:
    def __init__(self):
        self.file_info = {}  # Stores file details: {file_path: {'type': 'txt'/'fasta', 'lines'/ 'sequences': count, 'folder': folder_name}}
        self.lr_counts = {}  # Stores L/R counts for files: {file_path: {'Lxx': count, 'Rxx': count}}
    def count_L_R(self, file_path):
        """
        Reads a .txt file to count lines and analyze L(number)/R(number) patterns.
        Parameters:
            file_path (str): Path to the .txt file.
        """
        line_count = 0
        lr_pattern_counts = {}
        folder_name = os.path.dirname(file_path)
        with open(file_path, 'r') as file:
            for line in file:
                line_count += 1
                first_column = line.split()[0] if line.split() else ""
                for part in first_column.split('_'):
                    if part.startswith(('L', 'R')) and part[1:].isdigit():
                        lr_pattern_counts[part] = lr_pattern_counts.get(part, 0) + 1

        self.file_info[file_path] = {'type': 'txt', 'lines': line_count, 'folder': folder_name}
        self.lr_counts[file_path] = lr_pattern_counts
    def count_sequences_and_l_r_patterns(self, file_path):
        """
        Reads a .fasta file to count the number of sequences and analyze L(number)/R(number) patterns in seqID.
        Parameters:
            file_path (str): Path to the .fasta file.
        """
        sequence_count = 0
        lr_pattern_counts = {}
        folder_name = os.path.dirname(file_path)
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    sequence_count += 1
                    seqID = line[1:].strip()  # Extract seqID, removing '>' and whitespace
                    for part in seqID.split('_'):
                        if part.startswith(('L', 'R')) and part[1:].isdigit():
                            lr_pattern_counts[part] = lr_pattern_counts.get(part, 0) + 1

        self.file_info[file_path] = {'type': 'fasta', 'sequences': sequence_count, 'folder': folder_name}
        self.lr_counts[file_path] = lr_pattern_counts
