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
ls *fastq.gz | awk 'BEGIN{FS="_";OFS="_"}{print $1, $2}' | sort | uniq > TechReplicates

# create file with list of biological replicates
ls *fastq.gz | awk 'BEGIN{FS="_";OFS="_"}{print $2}' | sort | uniq > BioReplicates
```

### SortMeRNA

We'll pull PE and SE aligned / non-aligned reads.

#### Get unaligned reads to assemble mRNA

For each technical replicate, the following files will be output in the folder `${wd}/sortmerna`:

- `${sample}_sortmerna_aligned_mRNA.allreads.PE.fastq`. Indeed these files are interleaved so maybe they need to be processed (de-interleaved) before running Trinity. If this is the case, you can load SortMeRNA and use the script `unmerge-paired-reads.sh`.
- `${sample}_sortmerna_aligned_mRNA.allreads.SE.fastq`.

**Run SortMeRNA**

- PE reads

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_getUnal.allreads.PE.%a \
-o sortmerna_getUnal.allreads.PE.%a.out \
-e sortmerna_getUnal.allreads.PE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# interleave fastq as required by SortMeRNA
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz > ${sample}.R1.fastq
zcat ${wd}/${sample}_f.aT.notCombined_2.fastq.gz > ${sample}.R2.fastq

merge-paired-reads.sh \
${sample}.R1.fastq \
${sample}.R2.fastq \
${sample}.interleaved.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5.8s-database-id98:$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5s-database-id98:$SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DBS/index/silva-arc-16s-id95:$SORTMERNA_DBS/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMERNA_DBS/index/silva-arc-23s-id98:$SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90:$SORTMERNA_DBS/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMERNA_DBS/index/silva-bac-23s-id98:$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95:$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DBS/index/silva-euk-28s-id98 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_mRNA.PE \
--other ${sample}_sortmerna_aligned_mRNA.allreads.PE \
--paired_in --paired_out --fastx --log \
--num_alignments 1 \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_mRNA.allreads.PE*
cp ${sample}_sortmerna_aligned_mRNA.allreads.PE* ${wd}/sortmerna
cp ${sample}_sortmerna_mRNA.PE*log ${wd}/sortmerna

BWE
```

- SE reads

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_getUnal.allreads.SE.%a \
-o sortmerna_getUnal.allreads.SE.%a.out \
-e sortmerna_getUnal.allreads.SE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# copy input file to tmp dir
zcat ${wd}/${sample}_f.aT.extendedFrags.fastq.gz > ${sample}.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5.8s-database-id98:$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5s-database-id98:$SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DBS/index/silva-arc-16s-id95:$SORTMERNA_DBS/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMERNA_DBS/index/silva-arc-23s-id98:$SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90:$SORTMERNA_DBS/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMERNA_DBS/index/silva-bac-23s-id98:$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95:$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DBS/index/silva-euk-28s-id98 \
--reads ${sample}.fastq \
--aligned ${sample}_sortmerna_mRNA.SE \
--other ${sample}_sortmerna_aligned_mRNA.allreads.SE \
--fastx --log \
--num_alignments 1 \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_mRNA.allreads.SE*
cp ${sample}_sortmerna_aligned_mRNA.allreads.SE* ${wd}/sortmerna
cp ${sample}_sortmerna_mRNA.SE*log ${wd}/sortmerna

BWE
```

#### Get read datasets for each domain/SSU (get aligned reads)

Need to interleave PE fastq files first!

- Bacterial SSU (PE reads)

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_bac.allreads.PE.%a \
-o sortmerna_bac.allreads.PE.%a.out \
-e sortmerna_bac.allreads.PE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# interleave fastq as required by SortMeRNA
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz > ${sample}.R1.fastq
zcat ${wd}/${sample}_f.aT.notCombined_2.fastq.gz > ${sample}.R2.fastq

merge-paired-reads.sh \
${sample}.R1.fastq \
${sample}.R2.fastq \
${sample}.interleaved.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_aligned_bacSSU.allreads.PE \
--paired_in --fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_bacSSU.allreads.PE*
cp ${sample}_sortmerna_aligned_bacSSU.allreads.PE* ${wd}/sortmerna
BWE
```

- Archaeal SSU (PE reads)

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_arc.allreads.PE.%a \
-o sortmerna_arc.allreads.PE.%a.out \
-e sortmerna_arc.allreads.PE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# interleave fastq as required by SortMeRNA
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz > ${sample}.R1.fastq
zcat ${wd}/${sample}_f.aT.notCombined_2.fastq.gz > ${sample}.R2.fastq

merge-paired-reads.sh \
${sample}.R1.fastq \
${sample}.R2.fastq \
${sample}.interleaved.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DBS/index/silva-arc-16s-id95 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_aligned_arcSSU.allreads.PE \
--paired_in --fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_arcSSU.allreads.PE*
cp ${sample}_sortmerna_aligned_arcSSU.allreads.PE* ${wd}/sortmerna
BWE
```

- Eukaryotic SSU (PE reads)

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_euk.allreads.PE.%a \
-o sortmerna_euk.allreads.PE.%a.out \
-e sortmerna_euk.allreads.PE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# interleave fastq as required by SortMeRNA
zcat ${wd}/${sample}_f.aT.notCombined_1.fastq.gz > ${sample}.R1.fastq
zcat ${wd}/${sample}_f.aT.notCombined_2.fastq.gz > ${sample}.R2.fastq

merge-paired-reads.sh \
${sample}.R1.fastq \
${sample}.R2.fastq \
${sample}.interleaved.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_aligned_eukSSU.allreads.PE \
--paired_in --fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_eukSSU.allreads.PE*
cp ${sample}_sortmerna_aligned_eukSSU.allreads.PE* ${wd}/sortmerna
BWE
```

- Bacterial SSU (SE reads)

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_bac.allreads.SE \
-o sortmerna_bac.allreads.SE.%a.out \
-e sortmerna_bac.allreads.SE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# copy files on tmp directory
zcat ${wd}/${sample}_f.aT.extendedFrags.fastq.gz > ${sample}.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90 \
--reads ${sample}.fastq \
--aligned ${sample}_sortmerna_aligned_bacSSU.allreads.SE \
--fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_bacSSU.allreads.SE*
cp ${sample}_sortmerna_aligned_bacSSU.allreads.SE* ${wd}/sortmerna
BWE
```

- Archaeal SSU (SE reads)

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_arc.allreads.SE \
-o sortmerna_arc.allreads.SE.%a.out \
-e sortmerna_arc.allreads.SE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# copy files on tmp directory
zcat ${wd}/${sample}_f.aT.extendedFrags.fastq.gz > ${sample}.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DBS/index/silva-arc-16s-id95 \
--reads ${sample}.fastq \
--aligned ${sample}_sortmerna_aligned_arcSSU.allreads.SE \
--fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_arcSSU.allreads.SE*
cp ${sample}_sortmerna_aligned_arcSSU.allreads.SE* ${wd}/sortmerna
BWE
```

- Eukaryotic SSU (SE reads)

```bash
cd ${wd}
mkdir -p sortmerna

sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < TechReplicates) \
-J sortmerna_euk.allreads.SE \
-o sortmerna_euk.allreads.SE.%a.out \
-e sortmerna_euk.allreads.SE.%a.err<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p TechReplicates)

# move to temporary directory
cd ${SNIC_TMP}

# copy files on tmp directory
zcat ${wd}/${sample}_f.aT.extendedFrags.fastq.gz > ${sample}.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95 \
--reads ${sample}.fastq \
--aligned ${sample}_sortmerna_aligned_eukSSU.allreads.SE \
--fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_eukSSU.allreads.SE*
cp ${sample}_sortmerna_aligned_eukSSU.allreads.SE* ${wd}/sortmerna
BWE
```
