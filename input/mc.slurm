#!/bin/bash
#SBATCH --job-name=diamond
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=4:00:00
#SBATCH --export=ALL

module purge
ml iccifort/2020.1.217 impi/2019.7.217 imkl/2020.1.217 

export FI_PROVIDER=tcp
srun vmc++.exe neel.json neel > out.txt 
