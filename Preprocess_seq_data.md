#### Working directory
`/proj/b2016308/nobackup/projects/viraja`.

### Raw reads QC (fastqc)

```bash
mkdir -p /pica/v8/b2013127_nobackup/projects/domeni/aspo/RNA2/fastqc
cd /pica/v8/b2013127_nobackup/projects/domeni/aspo/RNA2/fastqc
sbatch fastqc_wrapper.serial.sh \
"$(ls /proj/b2013127/nobackup/projects/domeni/aspo/RNA2/raw_reads/*fastq.gz)" .
```
with script `fastqc_wrapper.serial.sh` being:

```bash
#!/bin/bash

# should be used as
# fastqc_wrapper_serial.sh regexp_for_files outdir, eg
# sbatch fastqc_wrapper.serial.sh "$(ls /proj/b2013127/INBOX/M.Dopson_16_01/*/*/*.fastq.gz)" /proj/b2013127/nobackup/projects/domeni/aspo/RNA/

#SBATCH -A b2013127
#SBATCH -J fastqc
#SBATCH -o fastqc.out
#SBATCH -e fastqc.err
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 40:00
#SBATCH --mail-type=ALL
#SBATCH --mail-type=domenico.simone@lnu.se

module load bioinfo-tools
module load FastQC/0.10.1

infile_regex=$1
outdir=$2

time fastqc -t $(nproc) -o ${outdir} ${infile_regex}
```

### Raw reads low quality end and adapter trimming (trimmomatic)
```bash
cd /pica/v8/b2013127_nobackup/projects/domeni/aspo/RNA/raw_reads
ls *1.fastq.gz > dataset_file_trimmomatic
sbatch --array=1-$(wc -l < dataset_file_trimmomatic) --time=40:00 trimmomatic.sh
```

with script `trimmomatic.sh` being:

```bash
#!/bin/bash

# put in command line: array, time

#SBATCH -A b2013127 -p node
#SBATCH -J trimmomatic_%a -o trimmomatic_%a.out -e trimmomatic_%a.err
#SBATCH --mail-type=ALL --mail-user=domenico.simone@lnu.se

module load bioinfo-tools
module load trimmomatic/0.32

# eg
# infile: 1_101_1.fastq.gz
# infile_base: 1_101_filtered.fastq.gz
infile=$(sed -n "$SLURM_ARRAY_TASK_ID"p dataset_file_trimmomatic)
infile_base=$(echo ${infile} | awk 'BEGIN{FS="_";OFS="_"}{print $1, $2, "filtered.adapterTrimmed.fastq.gz"}')

time java -jar $TRIMMOMATIC_HOME/trimmomatic.jar PE \
-phred33 -threads $(nproc) -trimlog ${infile}.trimlog \
-basein ${infile} -baseout ${infile_base} \
ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:100
```
### Stats on filtering

`#grep ^Input trimmomatic_*.err | awk '{print $8, $13, $18, $21}'`
