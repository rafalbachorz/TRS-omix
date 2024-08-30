#!/bin/bash
#SBATCH --job-name=TRS_job            # Job name
#SBATCH --output=TRS_job_%j.out       # Standard output log
#SBATCH --error=TRS_job_%j.err        # Standard error log
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=16                # Number of CPU cores per task
#SBATCH --mem=32G                        # Memory per node
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec

source ~/miniconda3/etc/profile.d/conda.sh && conda activate TRS

# Path to your Python script
python_script="combined.py"

# Check if --help is passed as an argument
if [[ "$1" == "--help" ]]; then
    python $python_script --help
    exit 0
fi

# Run the Python script with all passed arguments
python $python_script "$@"
