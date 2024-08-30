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
4. **`TRS_part.py`**: This script is used to obtain initial results for subsequent BLAST analysis.

5. **`TRS_BLAST.sh`** BLASTS sequences from a specific directory against nt database at 100% identity threshold.

6. **`Blast_part.py`**: This script is used to obtain final results from BLAST files. 

7. **Combined pipeline using `combined.py`**: This script executes all the steps above including the BLASTing step and automatically detects which version of the script to use depending on slurm availability (as of now threads and memory parameters can be changed only in the scripts themselves)






6. Another column is added to each entry in the blast files containing acession matching TaxID
7. Ciąg dalszy .....

