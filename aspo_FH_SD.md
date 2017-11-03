### Setup working environment

Run this chunk every time you start or resume an analysis in this doc)

```bash
export wd="/pica/v9/b2016308_nobackup/projects/aspo_FH_SD"
export ordir="/pica/v8/b2013127_nobackup/projects/domeni/aspo"
```

Create symlinks of datasets and lists of BioReplicates / TechReplicates

```bash
cd $wd
ln -s ${ordir}/RNA*/processed_reads/*_?A????A-*_f.aT.extendedFrags.fastq.gz .
ln -s ${ordir}/RNA*/processed_reads/*_?A????A-*_f.aT.notCombined_?.fastq.gz .

# delete files for SA1229A-1R
for i in $(ls *SA1229A-1R*); do unlink $i; done

# create file with list of technical replicates
ls | awk 'BEGIN{FS="_";OFS="_"}{print $1, $2}' | sort | uniq > TechReplicates

# create file with list of biological replicates
ls | awk 'BEGIN{FS="_";OFS="_"}{print $2}' | sort | uniq > BioReplicates
```

### SortMeRNA

Run on all reference db simultaneously to get the **unaligned** output (reads used to do the mRNA coassembly)

#### Get read datasets for each domain/SSU

Need to interleave fastq files first!

- Bacterial SSU (PE reads)

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_${sample}_bac.allreads \
-o sortmerna_${sample}_bac.allreads.out \
-e sortmerna_${sample}_bac.allreads.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# interleave fastq as required by SortMeRNA
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz > ${sample}.R1.fastq
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz  > ${sample}.R2.fastq

merge-paired-reads.sh \
${sample}.R1.fastq \
${sample}.R2.fastq \
${sample}.interleaved.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_aligned_bacSSU.allreads \
--paired_in --fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_bacSSU.allreads*
cp ${sample}_sortmerna_aligned_bacSSU.allreads* ${wd}/sortmerna
BWE
```

- Archaea SSU

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_${sample}_arc.allreads \
-o sortmerna_${sample}_arc.allreads.out \
-e sortmerna_${sample}_arc.allreads.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# interleave fastq as required by SortMeRNA
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz > ${sample}.R1.fastq
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz  > ${sample}.R2.fastq

merge-paired-reads.sh \
${sample}.R1.fastq \
${sample}.R2.fastq \
${sample}.interleaved.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DBS/index/silva-arc-16s-id95 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_aligned_arcSSU.allreads \
--paired_in --fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_arcSSU.allreads*
cp ${sample}_sortmerna_aligned_arcSSU.allreads* ${wd}/sortmerna
BWE
```
