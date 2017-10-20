### Run Kaiju

```bash
#salloc -p devel -t 1:00:00 -A b2016308

cd /home/domeni/projects_b2016308/aspo
ls *.fa > samples

sbatch -t 5:00:00 -A b2013127 -p node \
--array=1-$(wc -l < samples) \
-J kaiju -e kaiju_%a.err -o kaiju_%a.out<<'EOF'
#!/bin/bash

module load bioinfo-tools
module load Krona
export PATH=/home/domeni/downloads/kaiju/bin/:$PATH

infile=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples)

time kaiju \
-z 16 \
-t /proj/b2016308/nobackup/dl/nodes.dmp \
-f /proj/b2016308/nobackup/dl/kaiju_db_nr_euk.fmi \
-i ${infile} \
-o ${infile}.kaiju
EOF
```

Mean run time: ~ 7' (n. of seqs: 818 - 25586).

**Generate reports**

Don't need to run as batch script indeed.

```bash
export kaijuData="/proj/b2016308/nobackup/dl/"

sbatch -t 4:00:00 -A b2013127 -p core \
--array=1-$(wc -l < samples) \
-J kaiju_post -e kaiju_post_%a.err -o kaiju_post_%a.out<<'EOF'
#!/bin/bash

module load bioinfo-tools
export PATH=/home/domeni/downloads/kaiju/bin/:$PATH

infile=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples)

# generate reports
kaijuReport \
-t ${kaijuData}/nodes.dmp \
-n ${kaijuData}/names.dmp \
-i ${infile}.kaiju \
-r genus \
-o ${infile}.kaiju.out.summary

addTaxonNames \
-t ${kaijuData}/nodes.dmp \
-n ${kaijuData}/names.dmp \
-i ${infile}.kaiju \
-o ${infile}.kaiju-names.out \
-p

kaiju2krona \
-t ${kaijuData}/nodes.dmp \
-n ${kaijuData}/names.dmp \
-i ${infile}.kaiju \
-o ${infile}.kaiju.out.krona

ktImportText \
-o ${infile}.kaiju.out.html \
${infile}.kaiju.out.krona
EOF
```
