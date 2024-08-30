# TRS-omix: Search Engine

This code accompanies the paper:

- **Sebastian Sakowski, Marta Majchrzak, Jacek Waldmajer, Pawel Parniewski**: *TRS-omix: a new search engine for trinucleotide flanked sequences*. 2021.

The content of the repository is derived from its predecessor, available at:
[https://github.com/TRS-omix/software](https://github.com/TRS-omix/software)

## ToDo's:

### Simple Tasks

1. **Sequence Flanking Correction**: Add flanking sequences (`*CGACGACGACG*`) analogously on the right side. 

### Advanced

1. **Introduce a way to resume processing from the last completed step**: Find crucial points in pipeline, after their completion add currently stored variables and info about present files to .json (or other format). Bonus points with args we should instantly know what the name of the folder *should be* so we can instantly do a search (function for it is present) and prompt the user for folder if we find multiple (we are using fnmatch) if .json is found load it. The problem here is that i have no experience with something like this so I'll need help with creating the logic behind it. 

# Operating Mechanism

> [!IMPORTANT]
> 1. **Environment Setup with Conda**: The necessary package list for script operation is in `environment.yml`.
>    
>    Quick installation command in terminal: `conda create env -f environment.yml -n TRS`

> [!IMPORTANT]
> 2. **Activate Environment**: Use `conda activate TRS`. **All further operations should be performed in this environment**.

> [!WARNING]
> 3. **Compilation of TRS-wrapper**: (Request to Mr. Rafal for the exact script needed for compilation)

> [!NOTE]
> 4. **Usage of `TRS_part.py`**: This script is used to obtain initial results for subsequent BLAST analysis. Detailed operation described below.Or click [HERE](#trs-and-fasta)
> 5. **Proceed with BLASTING the obtained** `.fasta` sequences against nt database with tabular output format and 100% identity. Using TRS_BLAST.sh (slurm version)
> 6. **Usage of `Blast_part.py`**: This script is used to obtain final results from blast files
> 7. **Combined pipeline using `combined.py`**: This script executes all the steps above including the BLASTing step and automatically detects which version of the script to use depending on slurm availability (as of now threads and memory parameters can be changed only in the scripts themselves)

> [!CAUTION]
> 8. **DO NOT REMOVE/MOVE THE DIRECTORY CREATED AFTER RUNNING THE SCRIPTS**

## TRS and fasta
usage: TRS_part.py [-h] --input_fasta_folder_path INPUT_FASTA_FOLDER_PATH
                   --tmin TMIN --tmax TMAX --mode MODE [--redo] [--cont]
                   --email EMAIL --threshold THRESHOLD --length_to_extract
                   LENGTH_TO_EXTRACT [--cd_hit_path CD_HIT_PATH]

This program extracts TRS sequences from a series of input genomes, allows for
length selection of extracted fragments, and prepares sequences for further
analysis.

optional arguments:
  -h, --help            show this help message and exit
  --input_fasta_folder_path INPUT_FASTA_FOLDER_PATH
                        Path to a folder containing genomes in fasta format
                        from which TRS sequences will be extracted[REQUIRED]
  --tmin TMIN           Minimum length of TRS sequences[REQUIRED]
  --tmax TMAX           Maximum length of TRS sequences[REQUIRED]
  --mode MODE           Mode of operation, must be 0 or 1[REQUIRED]
  --redo                Redo the analysis if results directory already exists[NOT IMPLEMENTED]
  --cont                Continue the analysis from saved TRS results file[NOT IMPLEMENTED]
  --email EMAIL         Address e-mail to be used for connection with NCBI[REQUIRED]
                        databases
  --threshold THRESHOLD
                        Identity threshold for clustering using cdhit has to
                        be between 0.8 and 1.0[REQUIRED]
  --length_to_extract LENGTH_TO_EXTRACT
                        Length of flanking sequences to be extracted from the
                        full TRS sequence[REQUIRED]
  --cd_hit_path CD_HIT_PATH
                        Path to the cd-hit-est executable
> [!IMPORTANT]
> 3. Creates a new path in the folder where the script is located, named according to the pattern:
>    
>    `inputdirectory_results`
>    
>    This path will contain all files generated during analysis, dynamically changing to include information about the experiment.

> [!CAUTION]
> 4. **Do not modify this folder in _any_ way**. Files within it can be copied elsewhere, but the original location and file names must be preserved—at least for now.

> [!NOTE]
> 5. Upon specifying required parameters and creating the folder, TRS-omix operates, typically taking 1-1.5 hours for 7 genomes. The analysis duration also depends on the maximum length parameter.

> [!IMPORTANT]
> 6. TRS output is saved in the `TRS_output` folder, with an equivalent to `interiors.txt` from TRS-omix csv file named:
>    
>    `inputdirectory_results.csv`

7. Users are asked for their email, which will be used to obtain organism names along with the `GENOME` column of `inputdirectory_results.csv`. It's best if the `.fasta` files originate from NCBI Nucleotide for compatibility.

8. Extracted sequence fragments are saved to a `.fasta` file named "combined_sequences.fasta", with sequences named according to the scheme: `Species_name_L/R{number}`. The number accompanying L/R are encoded trinucleotide repeats.

9. Each L/R pair receives a number indicating a pair `Species_name_L/R{number}_{pair_number}`, and sequences are saved to `combined_sequences_unique.fasta`.

> [!IMPORTANT]
> 10. Subsequent operations include clustering with cd-hit (automated if cd-hit is installed or will prompt for path if not found) and setting the desired identity degree.
>     Note that this is one of the most time-consuming processes in the current script but is highly dependent on the desired identity threshold, being longest for 0.75 and **_very_ short** for 1.0.

11. The script also performs operations on clusters to clean them and obtain sequence IDs to be discarded. cd-hit results are located in the `cd-hit results` folder.

12. Two new fasta files are created in the `filtered_sequences` folder, one containing sequences within clusters and another outside them.

> [!WARNING]
> 13. These files should then be BLASTed against the nt database with parameters `perc_identity 100` and `-outfmt 6`.


## BLAST processing

> [!CAUTION]
> 1. Remember to put your blast output files into blast_output directory!

2. Script will attempt to find the blast_output folder on the users machine if multiple are found user will be prompted for selection

> [!NOTE]
> 3. Extensions will be added to all present files(so keep only blast files here) and unique acession - seqID pairs will be found

4. Further processing will require user to obtain *nucl_gb.accession2taxid* database file from NCBI FTP available at : https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/ download and extract the file. 

5. This database will help us get a quick list of all acession - taxid pairs in our dataset

6. Another column is added to each entry in the blast files containing acession matching TaxID
7. Ciąg dalszy .....

