# Run Interproscan on (meta)genome assembly

## Extract ORFs with prokka

Use script in this repo, `run_prokka.sh`.

```bash
export inputFile="thawponds_assembly.fa"
export outName="THAWP"
mkdir -p logs

sbatch -A snic2018-3-22 -p node -t 20:00:00 \
-J prokka -o logs/prokka_thawponds.out -e logs/prokka_thawponds.err \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

# SBATCH -J prokka_scilife -o prokka_scilife_%a.out -e prokka_scilife_%a.err
# SBATCH --array=1-240
# SBATCH -p node -t 2:30:00 -A b2016308
# SBATCH --mail-user=ALL

module load prokka

#inputFile=$(sed -n "$SLURM_ARRAY_TASK_ID"p binList)
#outName=$(echo ${inputFile} | sed 's/.fa//g') 

echo $inputFile
echo $outName

prokka --outdir ${outName} \
--locustag ${outName} \
--cpus 16 \
--norrna \
--notrna \
--evalue 1e-20 \
${inputFile}

EOF
```

