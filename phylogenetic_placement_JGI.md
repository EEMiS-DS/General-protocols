# Phylogenetic placement of NGS reads

### PaPaRa

#### Workflow from SSU read sorting to phylogenetic placement

**Dataset**: JGI CSP metagenome data.

Tools needed:
- SortMeRNA (available on cluster)
- PaPaRa (in `/proj/b2016308/glob/` directory: the precompiled executable has problems, so had to re-compile it)
- RAxML (available on cluster)

Files needed:
- one phylogenetic tree for each domain of life
- MSAs from which the phylogenetic trees were computed

#### Define environment variables

```bash
# folder with shared tools, scripts etc.
export PATH=/proj/b2016308/glob:$PATH

# working directory
export wd="/pica/v9/b2016308_nobackup/projects/JGI_CSP_analyses/phylogenetic_placement"
export samples="11383.1.204469.GTGAAA 11383.1.204469.GTGGCC 11408.1.205223.GTTTCG 11341.6.202084.TAGCTT 11606.7.214304.GTCCGC 11287.7.199536.GGCTAC 11287.8.199539.CTTGTA 11287.8.199539.AGTCAA 11292.4.199689.AGTTCC 11292.4.199689.ATGTCA 11292.5.199692.CCGTCC 11292.5.199692.GTAGAG"
# reference alignment in phylip format
export arcSSU_RA="/home/domeni/projects_b2016308/TOL/170921/TOS_all.l600.ark.clean.95Gaps.afa.reduced"
export bacSSU_RA="/home/domeni/projects_b2016308/TOL/170921/TOS_all.l600.bac.clean.95Gaps.afa.reduced"
# reference tree in newick format
export arcSSU_RT="/home/domeni/projects_b2016308/TOL/170921/RAxML_bipartitionsBranchLabels.TOS_all.l600.ark.clean.95Gaps.reduced_n4"
export bacSSU_RT="/home/domeni/projects_b2016308/TOL/170921/RAxML_bipartitionsBranchLabels.TOS_all.l600.bac.clean.95Gaps.reduced_n4"

# output folders
export sortmernaChunkFolder=${wd}/sortmerna_out_chunks
export paparaOutFolder=${wd}/papara_out_chunks
export raxmlEPAChunkFolder=${wd}/raxmlEPA_out_chunks
```

### Run analysis

#### Setup working environment

```bash
# create folder
mkdir -p ${wd}
cd ${wd}
ln -s /proj/b2016308/nobackup/projects/JGI_CSP_data/*METAGENOME* .
```

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
cd ${wd}
mkdir -p ${paparaOutFolder}

for domain in arc bac; do
#for domain in bac; do
    export domain=${domain}
    if [[ ${domain} == "arc" ]]
    then
        export RT=${arcSSU_RT}
        export RA=${arcSSU_RA}
        ls ${sortmernaChunkFolder}/*${domain}SSU*.fa | awk 'BEGIN{FS="/"}{print $NF}' > sortmerna_out_chunks.${domain}.list
    else
        export RT=${bacSSU_RT}
        export RA=${bacSSU_RA}
        ls ${sortmernaChunkFolder}/*${domain}SSU*.fa | awk 'BEGIN{FS="/"}{print $NF}' > sortmerna_out_chunks.${domain}.list
    fi
sbatch -t 10:00:00 -p node -A b2016308 \
--array=1-$(wc -l < sortmerna_out_chunks.${domain}.list) \
-J papara_${domain}_%a \
-o ${paparaOutFolder}/papara_${domain}_%a.out \
-e ${paparaOutFolder}/papara_${domain}_%a.err \
--mail-type=ALL --mail-user=domenico.simone@lnu.se<<'BWE'
#!/bin/bash

# test
# infile="P1607_145_sortmerna_aligned_bacSSU.allreads.orphans.1.fa"
infile=$(sed -n "$SLURM_ARRAY_TASK_ID"p sortmerna_out_chunks.${domain}.list)

export PATH=/proj/b2016308/glob/:$PATH
##for testing
#export sample=P1607_145
#export i=R1

cp ${RT} ${SNIC_TMP}
cp ${RA} ${SNIC_TMP}
cp ${sortmernaChunkFolder}/${infile} ${SNIC_TMP}
cd ${SNIC_TMP}

echo "temp folder:"
ls *

time papara \
-j 16 \
-t $(basename ${RT}) \
-s $(basename ${RA}) \
-q ${infile} \
-n ${infile}

tar -cvzf \
${infile}.papara.tar.gz *

cp \
${infile}.papara.tar.gz \
${paparaOutFolder}
BWE
done
```

Runtime for 100K seqs: ~ 6'30".

#### RAxML-EPA

```bash
mkdir -p ${raxmlEPAChunkFolder}

for domain in arc bac; do
    export domain=${domain}
    if [[ ${domain} == "arc" ]]
    then
        export RT=${arcSSU_RT}
        export depjid="11289611"
    else
        export RT=${bacSSU_RT}
        export depjid="11289612"
    fi
# use variables:
#   arcSSU_RT
#   sortmernaChunkFolder
#   paparaOutFolder which contains:
#       - archive ${infile}.papara.tar.gz with alignment file
#           - papara_alignment.${infile}
# for infile in $(ls $sortmernaChunkFolder | grep -v P1607_145_sortmerna_aligned_arcSSU.allreads.R1.0.fa); do
#     export infile=$infile
# if [[ ${infile} == *"bacSSU"* ]]; then
#     export depjid="11289612"
# elif [[ ${infile} == *"arcSSU"* ]]; then
#     export depjid="11289612"
# fi
#export infile="P1607_145_sortmerna_aligned_arcSSU.allreads.R1.0.fa"

#salloc -p devcore -n 8 -t 1:00:00 -A b2013127

sbatch -p core -n 8 -t 20:00:00 -A b2016308 \
--array=1-$(wc -l < sortmerna_out_chunks.${domain}.list) \
--dependency=afterany:${depjid} \
-J psy_raxmlEPA_${domain}_%a \
-o ${raxmlEPAChunkFolder}/raxmlEPA_${domain}_%a.out \
-e ${raxmlEPAChunkFolder}/raxmlEPA_${domain}_%a.err \
--mail-type=ALL --mail-user=dome.simone@gmail.com<<'EOF'
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

#### Run all on total datasets



Archaea from
