#!/bin/bash

#SBATCH -J prokka_scilife -o prokka_scilife_%a.out -e prokka_scilife_%a.err
# SBATCH --array=1-240
#SBATCH -p node -t 2:30:00 -A b2016308
#SBATCH --mail-user=ALL --mail-user=annalisa.16@hotmail.it

module load prokka

inputFile=$(sed -n "$SLURM_ARRAY_TASK_ID"p binList)
outName=$(echo ${inputFile} | sed 's/.fa//g') 

echo $inputFile
echo $outName

prokka --outdir ${outName} \
--locustag ${outName} \
--cpus 16 \
--norrna \
--notrna \
--evalue 1e-20 \
${inputFile}
