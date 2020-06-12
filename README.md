# wastewater_sarscov2
A repository for sharing data, analyses, and results from Bay Area wastewater Sars-CoV-2 sequencing

Internal note: The location for this repository and additional data on the IGI compute cluster is:
`/groups/banfield/projects/industrial/wastewater_bayarea/2020`

The strategy for analyzing samples currently is three parted:

## 1. Metagenomic assembly with megahit

Example command:
`megahit -1 ./052920_01/5_13_A_S1_R1_001.fastq.gz -2 ./052920_01/5_13_A_S1_R2_001.fastq.gz -t 48 -o 052920_01_asm`

We can then call protein predictions in our assembly:
`prodigal -i 052920_01.contigs.fa -a 052920_01.contigs.faa -p meta`

I would then blast these against a database of, e.g. human viral proteins from NCBI:

```
blastp -query 052920_01.contigs.faa -db ../../viral_db/human_viruse_proteins.faa -outfmt "6 qacc pident evalue qlen sacc stitle" -evalue 0.0001 > assembly_blast.output
```

We can also do much larger BLASTs, including the BLAST to uniprot and KEGG, as is done for ggkbase import, later.

In the future with NovaSeq runs, we will do metagenomic binning and genome curation, of course.

## 2. Mapping to the SARS-CoV-2 reference

(`./reference.fna`)
`bowtie2 -p 20 -x reference.fasta -1 ../raw/052920_04/5_19_A_S4_R1_001.fastq.gz -2 ../raw/052920_04/5_19_A_S4_R2_001.fastq.gz  | sambam > 052920_04.bam`

`sambam` is a custom samtools script on the compute cluster, reproduced below:
```
#!/bin/sh

#v0.1: Simple wrapper to convert sam files to sorted bam files - Rohan Sachdeva

/shared/software/samtools/latest/bin/samtools view -bu | /shared/software/samtools/latest/bin/samtools sort -@6 -l 6 -
```

We can then remove PCR duplicates from this BAM (often not a big concern, but seriously important for these low input libraries):
```
./sambamba-0.7.1-linux-static markdup -r ../052920_04.bam 052920_04_nodup.bam
```

We then profile this BAM with inStrain:
```
inStrain profile 052920_04_nodup.bam -l 0.9 ../reference.fasta -o 052920_04_instrain
```

We will want to manually verify and filter every SNP. The next version of inStrain (hopefully out soon) should have some SNP filtering that is high quality enough to believe every SNP without this, but for this project we will still want to look at every SNP we find.

## 3. Known viral analysis with Kraken-uniq
Kraken-uniq paper:
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0

In general, these metagenomic classifiers have *extremely* poor specificity and I do not condone their use (no matter how many citations they have) in most circumstances. However, a properly filtered Kraken-uniq result is still sometmies useful. The key thing to filter on is the 'cov' column, which is breadth of sequencing. Anything with cov < 0.5 (50%) is near useless, filter to cov > 0.5 in the reports. Here we try to classify any known viral sequences:

```
../../krakenuniq/krakenuniq --db ../../krakenuniq/kraken_viral_db/ --threads 16 --report-file 052920_04.report ../raw/052920_04/5_19_A_S4_R1_001.fastq.gz ../raw/052920_04/5_19_A_S4_R2_001.fastq.gz
```

I still really wouldn't publish results from kraken-uniq, but it is a nice first pass to obtain some information about known viruses in your sample. Note that it identifies a bocavirus (in 052920_04.report) that is also in the assembly on that sample. 

