# Run hmmer (ref db = Pfam) on (meta)genome assembly

## Annotate with prokka

You don't really need it but it's ok.
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

## Get ORFs with prodigal

```bash
cd /home/domeni/thaw_ponds_contig_mappings/new_contig_mappings

mkdir -p logs
mkdir -p prodigal

sbatch -p core -t 10:00:00 -A snic2018-3-22 \
-J prodigal -o logs/prodigal_thaw_ponds.out -e logs/prodigal_thaw_ponds.err \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load prodigal

prodigal -i thawponds_assembly.fa \
-a prodigal/thawponds_assembly.cds.faa \
-d prodigal/thawponds_assembly.cds.ffn \
-o prodigal/thawponds_assembly.cds.out \
-p meta \
-m

EOF
```

## Test execution time on 100 sequences

```bash
mkdir test_hmmer
head -n200 THAWP/sprot.faa > test_hmmer/test.faa
cd test_hmmer

sbatch -A snic2018-3-22 -p node -t 20:00:00 \
-J hmmsearch_test -o hmmsearch_test.out -e hmmsearch_test.err \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load hmmer

time hmmsearch \
--tblout test.tblout \
-E 1e-5 \
--cpu 20 \
/home/domeni/thaw_ponds/pfam_db_2018/Pfam-A.hmm \
test.faa
EOF
```

Time for 100 sequences: 1'. We can split the dataset (n=2638865) in chunks of 60,000 sequences which should take 10 hours each.

## Split sequences and run hmmsearch

```bash
splitSeqFile.py <infile> \
fasta \
fasta \
60000

hmmsearch \
--tblout <output file> \
-E 1e-5 \
--cpu 20 \
~/Pfam/Pfam-mobility.hmm \
<input file (protein format)> > /dev/null
```
