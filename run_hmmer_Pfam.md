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
-f gff \
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

Time for 100 sequences: 1'. We can split the dataset (n=2638865) in chunks of 30,000 sequences which should take ~5 hours each.

## Split sequences and run hmmsearch

```bash
export wdir=`pwd`

splitSeqFile.py prodigal/thawponds_assembly.cds.faa \
fasta \
fasta \
30000

ls prodigal/thawponds_assembly.cds.*faa > prodigal/thawponds_assembly.cds.faa.files

sbatch -A snic2018-3-22 -p node -t 20:00:00 \
-J hmmer -o logs/hmmer_thawponds_%a.out -e logs/hmmer_thawponds_%a.err \
--array=1-$(wc -l < prodigal/thawponds_assembly.cds.faa.files) \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load hmmer

inputFile=$(sed -n "$SLURM_ARRAY_TASK_ID"p prodigal/thawponds_assembly.cds.faa.files)
basenameFile=$(basename $inputFile)
outFile=${basenameFile/.faa/.hmmer_pfam.tblout}
cp $inputFile ${SNIC_TMP}
cp /home/domeni/thaw_ponds/pfam_db_2018/Pfam-A.hmm ${SNIC_TMP}

cd ${SNIC_TMP}

hmmsearch \
--tblout $outFile \
-E 1e-5 \
--cpu 20 \
Pfam-A.hmm \
${basenameFile}

cp $outFile $wdir

EOF
```

Indeed 30000 sequences take less than 20' :-) How long does it take for the entire dataset (n=2638865) to run?

```bash
cd /home/domeni/thaw_ponds_contig_mappings/new_contig_mappings
mkdir -p hmmer
export wdir=`pwd`

sbatch -A snic2018-3-22 -p node -t 50:00:00 \
-J hmmer_all -o logs/hmmer_thawponds_all.out -e logs/hmmer_thawponds_all.err \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load hmmer

inputFile=prodigal/thawponds_assembly.cds.faa
basenameFile=$(basename $inputFile)
outFile=${basenameFile/.faa/.all.hmmer_pfam.tblout}
cp $inputFile ${SNIC_TMP}
cp /home/domeni/thaw_ponds/pfam_db_2018/Pfam-A.hmm ${SNIC_TMP}

cd ${SNIC_TMP}

hmmsearch \
--tblout $outFile \
-E 1e-5 \
--cpu 20 \
Pfam-A.hmm \
${basenameFile}

cp $outFile $wdir/hmmer

EOF
```
