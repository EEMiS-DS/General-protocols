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

## Perform read count for each predicted cds

Use `htseq-count` to perform two different read counts:

- one considering **all** (also **non-unique**) read placements. These should be input for things like **descriptive plots**;
- one considerin **only** unique read placements. These should be input for downstream analysis such as **differential gene expression** (.

### Test htseq-count for all counts

```bash
# test on smallest align file
#D-1_S32_L006.sorted.bam

# salloc -p devcore -t 1:00:00 -A snic2018-3-22

module load bioinfo-tools

module load htseq

time htseq-count \
-f bam \
-r name \
-s no \
-t CDS \
-i ID \
-m union \
-o D-1_S32_L006.counts.all.sam \
--nonunique all \
D-1_S32_L006.sorted.bam \
prodigal/thawponds_assembly.cds.out > counts/D-1_S32_L006.counts.all.out
```

### Test htseq-count for unique counts

```bash
# test on smallest align file
#D-1_S32_L006.sorted.bam

# salloc -p devcore -t 1:00:00 -A snic2018-3-22

module load bioinfo-tools
module load htseq

htseq-count \
-f bam \
-r pos \
-s no \
-t CDS \
-i ID \
-m union \
--nonunique none \
D-1_S32_L006.sorted.bam \
prodigal/thawponds_assembly.cds.out

```

### Run samtools sort + htseq-count (nonunique) + htseq-count (all)

With htseq-count, the best is to have bam sorted by read name. So we'll sort them by name and run htseq

```bash
export wdir=`pwd`

#mkdir -p alignments
ls alignments/*.sorted.bam | sed 's/.sorted.bam//g' | awk 'BEGIN{FS="/"}{print $2}' > sampleList

sbatch -A snic2018-3-22 -p core -t 20:00:00 \
-J ssort_htseq_thawponds -o logs/ssort_htseq_thawponds_%a.out -e logs/ssort_htseq_thawponds_thawponds_%a.err \
--array=1-$(wc -l < sampleList) \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load pysam/0.13-python2.7.11
module load htseq
module load samtools

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p sampleList)

cp prodigal/thawponds_assembly.cds.out ${SNIC_TMP}
cp alignments/${sample}.sorted.bam ${SNIC_TMP} && cd ${SNIC_TMP}

samtools sort -o ${sample}.sorted.name.bam ${sample}.sorted.bam

cp ${sample}.sorted.name.bam ${wdir}/alignments

time htseq-count \
-f bam \
-r name \
-s no \
-t CDS \
-i ID \
-m union \
--nonunique all \
${sample}.sorted.name.bam \
thawponds_assembly.cds.out > ${sample}.counts.all.out

time htseq-count \
-f bam \
-r name \
-s no \
-t CDS \
-i ID \
-m union \
--nonunique none \
${sample}.sorted.name.bam \
thawponds_assembly.cds.out > ${sample}.counts.unique.out

cp *counts*out ${wdir}/counts

EOF
```

Then concatenate count files, adding to each row the sample name and compressing the resulting file

```bash
catFile=all_samples.counts.all.out
rm -f $catFile
for i in $(ls *.all.out); do
    sampleName=$(echo $i | cut -f1 -d_)
    echo $sampleName
    awk -v sampleName="$sampleName" 'BEGIN{OFS="\t"}{print sampleName, $0}' $i >> $catFile 
done
gzip $catFile

catFile=all_samples.counts.unique.out
rm -f $catFile
for i in $(ls *.unique.out); do
    sampleName=$(echo $i | cut -f1 -d_)
    echo $sampleName
    awk -v sampleName="$sampleName" 'BEGIN{OFS="\t"}{print sampleName, $0}' $i >> $catFile 
done
gzip $catFile
```

## Test hmmer execution time on 100 sequences

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

## Split sequences (chunk size = 30000) and run hmmsearch

```bash
export wdir=`pwd`

splitSeqFile.py prodigal/thawponds_assembly.cds.faa \
fasta \
fasta \
30000

mkdir -p prodigal/split30000
mv prodigal/thawponds_assembly.cds.split30000.*.faa prodigal/split30000
ls prodigal/split30000/thawponds_assembly.cds.split30000.*.faa > prodigal/split30000/thawponds_assembly.cds.faa.files

export outDir=${wdir}/hmmer/split30000
mkdir -p ${outDir}

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

cp $outFile $outDir

EOF
```

Then concatenate results to be exported:
- hmmer results, filtering out comment lines and compressing the output;
- read counts with ambiguous reads filtered out (for differential abundance analysis);
- read counts with all reads (for descriptive plots).

```bash
# hmmer results
cd hmmer
grep -hv "^#" thawponds_assembly.cds.*.hmmer_pfam.tblout | gzip - > thawponds_assembly.cds.all.hmmer_pfam.tblout.gz
cd ..

# read counts with ambiguous reads filtered out

```

## Test hmmer on single file

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

Execution time: ~28h.

## Split sequences (chunk size = 500000) and run hmmsearch

What we need is the hmmer output got from the option `--domtblout`. But we will re-run hmmer on chunks of 500000 sequences.

```bash
export wdir=`pwd`

splitSeqFile.py prodigal/thawponds_assembly.cds.faa \
fasta \
fasta \
500000

mkdir -p prodigal/split500000
mv prodigal/thawponds_assembly.cds.split500000.*.faa prodigal/split500000
ls prodigal/split500000/thawponds_assembly.cds.split500000.*.faa > prodigal/split500000/thawponds_assembly.cds.faa.files

export outDir=${wdir}/hmmer/split500000
mkdir -p ${outDir}

sbatch -A snic2018-3-22 -p node -t 20:00:00 \
-J hmmer_%a -o logs/hmmer_thawponds_split500000_%a.out -e logs/hmmer_thawponds_split500000_%a.err \
--array=1-$(wc -l < prodigal/split500000/thawponds_assembly.cds.faa.files) \
--mail-type=ALL --mail-user=domenico.simone@slu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load hmmer

inputFile=$(sed -n "$SLURM_ARRAY_TASK_ID"p prodigal/split500000/thawponds_assembly.cds.faa.files)
basenameFile=$(basename $inputFile)
outFile=${basenameFile/.faa/.hmmer_pfam.out}
pfamOutFile=${basenameFile/.faa/.hmmer_pfam.pfamtblout}
domOutFile=${basenameFile/.faa/.hmmer_pfam.domtblout}
tblOutFile=${basenameFile/.faa/.hmmer_pfam.tblout}
cp $inputFile ${SNIC_TMP}
cp /home/domeni/thaw_ponds/pfam_db_2018/Pfam-A.hmm ${SNIC_TMP}

cd ${SNIC_TMP}

hmmsearch \
-o $outFile \
--pfamtblout $pfamOutFile \
--domtblout $domOutFile \
--tblout $tblOutFile \
-E 1e-5 \
--cpu 20 \
Pfam-A.hmm \
${basenameFile}

cp $outFile $outDir
cp $pfamOutFile $outDir
cp $domOutFile $outDir
cp $tblOutFile $outDir

EOF
```

Running time = ~3h30'.

Concatenate results and rsync to local machine.

```{bash}
cat hmmer/split500000/thawponds_assembly.cds.split500000.0000*domtblout > hmmer/split500000/thawponds_assembly.cds.split500000.all.hmmer_pfam.domtblout
```
