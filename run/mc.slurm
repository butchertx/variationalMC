#!/bin/bash
#SBATCH --job-name=su3_triangle
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=500M
#SBATCH --time=00:10:00
#SBATCH --export=ALL

module purge
module load iccifort/2020.1.217 impi/2019.7.217 imkl/2020.1.217 
srun triangle_csl.exe csl_pi3.json
