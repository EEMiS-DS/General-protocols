# Phylogenetic placement of NGS reads

#### Workflow from SSU read sorting to phylogenetic placement

Tools needed and run:
- SortMeRNA (available on cluster)
- PaPaRa (in `/proj/sllstore2017037/glob` directory: the precompiled executable has problems, so had to re-compile it)
- RAxML (available on cluster)

Files needed:
- one phylogenetic tree for each domain of life
- MSAs from which the phylogenetic trees were computed

#### Define environment variables

**Run this every time you run any step of the pipeline!!!**

Copy `pp_setup.sh` file in the analysis where the analysis will be performed, edit the `samples` variable and source it.
Indeed you need to edit the `samples` variable only if you are running SortMeRNA.

```bash
export PATH=/proj/sllstore2017037/glob:$PATH
cp /proj/sllstore2017037/glob/pp_setup.sh .

# edit the samples variable...

source pp_setup.sh
```

### Run analysis

PaPaRa needs the reads to be already oriented in the same strand they're going to be aligned. Reads are already reversed and complemented (if needed) in the sam output of SortMeRNA.

#### SortMeRNA

Don't need to interleave fastq files first!

- Bacteria

```bash
cd ${wd}
for sample in ${samples}; do
    export sample=${sample}
sbatch -t 10:00:00 -p node -A b2016308 \
-J sortmerna_${sample}_bac.allreads \
-o sortmerna_${sample}_bac.allreads.out \
-e sortmerna_${sample}_bac.allreads.err \
--mail-type=ALL --mail-user=tyuamail@mail.com,domenico.simone@lnu.se<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

# move to temporary directory
cd ${SNIC_TMP}

zcat ${wd}/${sample}.filter-METAGENOME.fastq.gz > ${sample}.interleaved.fastq

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DBS/index/silva-bac-16s-id90 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_aligned_bacSSU.allreads \
--paired_in --fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_bacSSU.allreads*
cp ${sample}_sortmerna_aligned_bacSSU.allreads* $wd
BWE
done
```

- Archaea

```bash
cd ${wd}
for sample in ${samples}; do
    export sample=${sample}
sbatch -t 10:00:00 -p node -A b2016308 \
-J sortmerna_${sample}_arc.allreads \
-o sortmerna_${sample}_arc.allreads.out \
-e sortmerna_${sample}_arc.allreads.err \
--mail-type=ALL --mail-user=domenico.simone@lnu.se<<'BWE'
#!/bin/bash

module load bioinfo-tools
module load SortMeRNA/2.1b

# move to temporary directory
cd ${SNIC_TMP}

zcat ${wd}/${sample}.filter-METAGENOME.fastq.gz > ${sample}.interleaved.fastq

sortmerna --ref /sw/apps/bioinfo/SortMeRNA/2.1b/milou/sortmerna/rRNA_databases/silva-arc-16s-id95.fasta,/sw/apps/bioinfo/SortMeRNA/2.1b/milou/sortmerna/index/silva-arc-16s-id95 \
--reads ${sample}.interleaved.fastq \
--aligned ${sample}_sortmerna_aligned_arcSSU.allreads \
--paired_in --fastx --log \
--num_alignments 1 \
--sam \
-a 16 -e 1e-20

ls ${sample}_sortmerna_aligned_arcSSU.allreads*
cp ${sample}_sortmerna_aligned_arcSSU.allreads* $wd
BWE
done
```

- Get statistics of read alignments

```bash
cd $wd
for sample in ${samples}; do
    for domain in arc bac; do
        totalAligned=$(grep -vc "^@" ${sample}_sortmerna_aligned_${domain}SSU.allreads.sam)
        orphanAl=$(cat ${sample}_sortmerna_aligned_${domain}SSU.allreads.sam | awk '{count[$1]++}END{ for (j in count) print j, count[j] }' | awk '$2!=2{print $0}' | wc -l)
        let "PE=($totalAligned-$orphanAl)/2"
        echo ${sample} ${domain}
        echo "Total reads aligned: $totalAligned"
        echo "Orphan reads: $orphanAl"
        echo "Fragments with both reads aligned: $PE"
        echo "Fragments with one read aligned: $orphanAl"
done
done
```

- Process sam output with script `processSortMeRNAsam.chunks.py` to get:
    - R1 file
    - R2 file
    - orphan read file

```bash
cd ${wd}
mkdir -p ${sortmernaChunkFolder}
for sample in ${samples}; do
    export sample=${sample}
    for domain in arc bac; do
        export domain=${domain}
sbatch -p core -t 1:00:00 -A b2013127 -J sortProc.${sample}.${domain}.chunks -o sortProc.${sample}.${domain}.chunks.out -e sortProc.${sample}.${domain}.chunks.err --mail-type=ALL --mail-user=domenico.simone@lnu.se<<'EOF'
#!/bin/bash

processSortMeRNAsam.chunks.py \
${sample}_sortmerna_aligned_${domain}SSU.allreads.sam \
${sortmernaChunkFolder}
EOF
done
done
```

this script will generate subfolders with chunks of fasta files (100K seqs/file).

#### PaPaRa

Needs a reference alignment in phylip format.

```bash
# Run on both archaea and bacteria
pp_papara.sh "arc bac" 

# Only bacteria
pp_papara.sh "bac" 
```

Runtime for 100K seqs (archaea, 172 seqs in the tree): ~ 6'30".
Runtime for 80K seqs (archaea, 1639 seqs in the tree): ~ 30'.

#### RAxML-EPA

```bash
mkdir -p ${raxmlEPAChunkFolder}

for domain in arc bac; do
    export domain=${domain}
    if [[ ${domain} == "arc" ]]
    then
        export RT=${arcSSU_RT}
    else
        export RT=${bacSSU_RT}
    fi
#salloc -p devcore -n 8 -t 1:00:00 -A b2013127
sbatch -p core -n 8 -t 20:00:00 -A b2016308 \
--array=1-$(wc -l < sortmerna_out_chunks.${domain}.list) \
-J JGI_raxmlEPA_${domain}_%a \
-o ${raxmlEPAChunkFolder}/raxmlEPA_${domain}_%a.out \
-e ${raxmlEPAChunkFolder}/raxmlEPA_${domain}_%a.err<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load raxml

infile=$(sed -n "$SLURM_ARRAY_TASK_ID"p sortmerna_out_chunks.${domain}.list)

export paparaOutfile="${paparaOutFolder}/${infile}.papara.tar.gz"
export alignFile="papara_alignment.${infile}"

# extract papara output to SNIC_TMP
tar -xvzf ${paparaOutfile} -C ${SNIC_TMP} ${alignFile}
cp ${RT} ${SNIC_TMP}

cd ${SNIC_TMP}
raxmlHPC-PTHREADS-AVX -f v \
-s ${alignFile} \
-G 0.1 \
-t ${RT} \
-m GTRCAT \
-n ${alignFile} \
-T 8

tar -cvzf ${infile}.raxmlEPA.tar.gz RAxML_*

cp ${infile}.raxmlEPA.tar.gz ${raxmlEPAChunkFolder}
EOF
done
```
