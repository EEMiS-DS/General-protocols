### Setup working environment

Run this chunk every time you start or resume an analysis in this doc)

```bash
export wd="/pica/v9/b2016308_nobackup/projects/aspo_FH_SD"
```

Create symlinks of datasets

```bash
export ordir="/pica/v8/b2013127_nobackup/projects/domeni/aspo"

cd $wd
ln -s ${ordir}/RNA*/processed_reads/*_?A????A-*_f.aT.extendedFrags.fastq.gz .
ln -s ${ordir}/RNA*/processed_reads/*_?A????A-*_f.aT.notCombined_?.fastq.gz .

# delete files for SA1229A-1R
for i in $(ls *SA1229A-1R*); do unlink $i; done
```
