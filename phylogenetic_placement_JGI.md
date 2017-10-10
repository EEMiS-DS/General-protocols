# Phylogenetic placement of NGS reads

### PaPaRa

#### Test workflow from SSU read sorting to phylogenetic placement

**Dataset**: Archaea reads from acidophilic MFC metatranscriptome datasets (Gaofeng). Two replicates: P1607_145 and P1607_153.

Tools needed:
- SortMeRNA (available on cluster)
- PaPaRa (in my `glob` directory: the precompiled executable has problems, so had to re-compile it)
- RAxML (available on cluster)

Files needed:
- one phylogenetic tree for each domain of life
- MSAs from which the phylogenetic trees were computed

#### Define environment variables

```bash
# working directory
export wd="/pica/v9/b2016308_nobackup/projects/JGI_CSP_analyses/phylogenetic_placement"
export samples="11383.1.204469.GTGAAA 1383.1.204469.GTGGCC 11408.1.205223.GTTTCG 11341.6.202084.TAGCTT 11606.7.214304.GTCCGC 11287.7.199536.GGCTAC 11287.8.199539.CTTGTA 11287.8.199539.AGTCAA 11292.4.199689.AGTTCC 11292.4.199689.ATGTCA 11292.5.199692.CCGTCC 11292.5.199692.GTAGAG"
# reference alignment in phylip format
export arcSSU_RA="/home/domeni/projects_b2016308/TOL/170921/TOS_all.l600.ark.clean.95Gaps.afa.reduced"
# reference tree in newick format
export arcSSU_RT="/home/domeni/projects_b2016308/TOL/170921/RAxML_bipartitionsBranchLabels.TOS_all.l600.ark.clean.95Gaps.reduced_n4"

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

# JGI datasets come already interleaved, but we need to reformat headers from
# @HISEQ09:341:CA3UPANXX:7:1213:10106:44664 1:N:0:GGCTAC
# to
# @HISEQ09:341:CA3UPANXX:7:1213:10106:44664/1
python -c """
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
outhandle = open(outfile, 'w')

for i in open(infile, 'r'):
	if i.startswith('@HISEQ'):
		header = i.split()[0] + '/' + i.split()[1][0]
		outhandle.write(header + '\n')
	else:
		outhandle.write(i)

outhandle.close()
""" ${wd}/${sample}.filter-METAGENOME.fastq.gz ${sample}.interleaved.fastq

# # interleave fastq as required by SortMeRNA
# merge-paired-reads.sh \
# ${wd}/${sample}.filtered.adapterTrimmed_1P.fastq \
# ${wd}/${sample}.filtered.adapterTrimmed_2P.fastq \
# ${sample}.interleaved.fastq

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

# JGI datasets come already interleaved, but we need to reformat headers from
# @HISEQ09:341:CA3UPANXX:7:1213:10106:44664 1:N:0:GGCTAC
# to
# @HISEQ09:341:CA3UPANXX:7:1213:10106:44664/1
python -c """
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
outhandle = open(outfile, 'w')

for i in open(infile, 'r'):
	if i.startswith('@HISEQ'):
		header = i.split()[0] + '/' + i.split()[1][0]
		outhandle.write(header + '\n')
	else:
		outhandle.write(i)

outhandle.close()
""" ${wd}/${sample}.filter-METAGENOME.fastq.gz ${sample}.interleaved.fastq

# # interleave fastq as required by SortMeRNA
# merge-paired-reads.sh \
# ${wd}/${sample}.filtered.adapterTrimmed_1P.fastq \
# ${wd}/${sample}.filtered.adapterTrimmed_2P.fastq \
# ${sample}.interleaved.fastq

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
    totalAligned=$(grep -vc "^@" ${sample}_sortmerna_aligned_arcSSU.allreads.sam)
    orphanAl=$(cat ${sample}_sortmerna_aligned_arcSSU.allreads.sam | awk '{count[$1]++}END{ for (j in count) print j, count[j] }' | awk '$2!=2{print $0}')
    let "PE=($totalAligned-$orphanAl)/2"
    echo ${sample}
    echo "Fragments with both reads aligned: $PE"
    echo "Fragments with one read aligned: $orphanAl"
done
```

- Process sam output with script `processSortMeRNAsam.chunks.py` to get:
    - R1 file
    - R2 file
    - orphan read file

```bash
mkdir sortmerna_out_chunks
for sample in ${samples}; do
    export sample=${sample}
sbatch -p core -t 1:00:00 -A b2013127 -J sortProc.${sample}.chunks -o sortProc.${sample}.chunks.out -e sortProc.${sample}.chunks.err --mail-type=ALL --mail-user=domenico.simone@lnu.se<<'EOF'
#!/bin/bash

processSortMeRNAsam.chunks.py \
${sample}_sortmerna_aligned_arcSSU.allreads.sam \
sortmerna_out_chunks
EOF
done
```

this script will generate subfolders with chunks of fasta files (100K seqs/file).

#### PaPaRa
Needs a reference alignment in phylip format.

```bash
cd ${wd}
export sortmernaChunkFolder=${wd}/sortmerna_out_chunks
export paparaOutFolder=${wd}/papara_out_chunks
for infile in $(ls $sortmernaChunkFolder); do
    export infile=$infile
sbatch -t 60:00:00 -p node -A b2016308 \
-J papara_${infile} \
-o ${paparaOutFolder}/papara_${infile}.out \
-e ${paparaOutFolder}/papara_${infile}.err \
--mail-type=ALL --mail-user=domenico.simone@lnu.se<<'BWE'
#!/bin/bash

# papara is in my glob directory
##for testing
#export sample=P1607_145
#export i=R1

cp ${arcSSU_RT} ${SNIC_TMP}
cp ${arcSSU_RA} ${SNIC_TMP}
cd ${SNIC_TMP}

papara \
-j 16 \
-t ${arcSSU_RT} \
-s ${arcSSU_RA} \
-q ${sortmernaChunk}/${infile} \
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
export raxmlEPAChunkFolder=${wd}/raxmlEPA_out_chunks
mkdir -p ${raxmlEPAChunkFolder}
# use variables:
#   arcSSU_RT
#   sortmernaChunkFolder
#   paparaOutFolder which contains:
#       - archive ${infile}.papara.tar.gz with alignment file
#           - papara_alignment.${infile}
for infile in $(ls $sortmernaChunkFolder | grep -v P1607_145_sortmerna_aligned_arcSSU.allreads.R1.0.fa); do
    export infile=$infile
#export infile="P1607_145_sortmerna_aligned_arcSSU.allreads.R1.0.fa"
export paparaOutfile="${paparaOutFolder}/${infile}.papara.tar.gz"
export alignFile="papara_alignment.${infile}"
#salloc -p devcore -n 8 -t 1:00:00 -A b2013127

sbatch -p core -n 8 -t 20:00:00 -A b2016308 \
-J raxmlEPA_${infile} \
-o ${raxmlEPAChunkFolder}/raxmlEPA_${infile}.out \
-e ${raxmlEPAChunkFolder}/raxmlEPA_${infile}.err \
--mail-type=ALL --mail-user=domenico.simone@lnu.se<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load raxml

# extract papara output to SNIC_TMP
tar -xvzf ${paparaOutfile} -C ${SNIC_TMP} ${alignFile}
cp ${arcSSU_RT} ${SNIC_TMP}
#algFile="papara_alignment.papara_test_revcomp_edited"
#treeFile="RAxML_bipartitionsBranchLabels.TOS_all.l600.ark.clean.95Gaps.reduced_n4"

cd ${SNIC_TMP}
raxmlHPC-PTHREADS-AVX -f v \
-s ${alignFile} \
-G 0.1 \
-t ${arcSSU_RT} \
-m GTRCAT \
-n ${alignFile} \
-T 8

tar -cvzf ${infile}.raxmlEPA.tar.gz RAxML_*

cp ${infile}.raxmlEPA.tar.gz ${raxmlEPAChunkFolder}
EOF
```

#### Run all on total datasets



Archaea from
