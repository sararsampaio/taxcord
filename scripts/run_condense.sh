#!/usr/bin/env bash
# Example SLURM job that condenses an annotated BLAST table.
# Adjust the SBATCH directives and paths for your cluster, then submit with:
#   sbatch scripts/run_condense.sh <annotated_blast.txt> <condensed_output.txt>

#SBATCH --job-name=taxcord_condense
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=your.email@example.com
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=normal
#SBATCH --output=taxcord_condense_%j.out
#SBATCH --error=taxcord_condense_%j.err

set -euo pipefail

INPUT="${1:?usage: run_condense.sh <annotated_blast> <condensed_output>}"
OUTPUT="${2:?usage: run_condense.sh <annotated_blast> <condensed_output>}"

taxcord condense "$INPUT" "$OUTPUT"
