# TRS-omix: Search Engine

This code accompanies the paper:

- **Sebastian Sakowski, Marta Majchrzak, Jacek Waldmajer, Pawel Parniewski**: *TRS-omix: a new search engine for trinucleotide flanked sequences*. 2021.

The content of the repository is derived from its predecessor, available at:
[https://github.com/TRS-omix/software](https://github.com/TRS-omix/software)

## ToDo's:

### Simple Tasks

1. **Sequence Flanking Correction**: Add flanking sequences (`*CGACGACGACG*`) analogously on the right side. 

2. **Batch File for Genome List**: Upon program startup, the batch file should contain a list of genomes with arbitrary genome file names. **KIND OF SOLVED** CURRENT IMPLEMENTATION ASKS FOR PATH AND READS .FASTA FILES LOCATED THERE 

3. **FRG NO Column in Output**: Currently manually generated and needs automation. It should contain a string from the fasta file header (`>..." "`) with a prefix-index. **SOLVED**

4. **Sequence Similarity in `interiors.txt` ("Interiors")**: Address similarity of sequences within `interiors.txt`. **SOLVED**

5. **Run script with args instead of user input**: Current implementation has low scalability.

### Advanced

1. **Introduce a way to resume processing from the last completed step**: Find crucial points in pipeline, after their completion add currently stored variables and info about present files to .json (or other format). Bonus points with args we should instantly know what the name of the folder *should be* so we can instantly do a search (function for it is present) and prompt the user for folder if we find multiple (we are using fnmatch) if .json is found load it. The problem here is that i have no experience with something like this so I'll need help with creating the logic behind it. 

2. **Testing**: the current script was run on limited number of samples from klebsiella, avium, ecoli and citrobacter genomes. We need to test it's capabilities especially after introduction of automated dictionary creation and update. I wrote a short script that downloads a specified number of genomes from a given genus I will include it here. 

3. **Multithreading**: good idea would be to figure out how to multithread currently present function and test how it goes then rewrite for args.

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
> 4. **Usage of `TRS_and_fasta_revised.py`**: This script is used to obtain initial results for subsequent BLAST analysis. Detailed operation described below.Or click [HERE](#trs-and-fasta-revised)

5. **Proceed with BLASTING the obtained** `.fasta` sequences against nt database with tabular output format and 100% identity.

> [!CAUTION]
> 6. **DO NOT REMOVE/MOVE THE DIRECTORY CREATED BY `TRS_and_fasta_revised.py`**

7. **Move BLAST results** to the blast_output directory (should already be created) in the directory `TRS_and_fasta_revised.py` created.

8. **Usage of `BLAST_part_revised.py`**: This script ...

## TRS_and_fasta_revised

1. Asks the user for the location of the folder containing `.fasta` files. Full path specification is recommended.

2. Queries additional parameters to be used in TRS-omix (minimum and maximum length, mode).

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



