# wastewater_sarscov2
A repository for sharing data, analyses, and results from Bay Area wastewater Sars-CoV-2 sequencing.

Note: our internal SNV calling uses a 0-based index (so the first position in the genome is 0), so that indexing is used in this repository; we convert to 1-based in supplementary tables.

Internal note: The location for this repository and additional data on the IGI compute cluster is:
`/groups/banfield/projects/industrial/wastewater_bayarea/2020`

## 1. List of files

### Notebooks
`./notebooks/CovidWastewater_1_makeTables_3_allTables_v2.ipynb` - Creates tables from clinical data downloaded from GISAID.

`./notebooks/Primary_analysis_and_plots.ipynb` - Reproduce most figures in the manuscript.

`./notebooks/CovidWastewater_2_plotting_6.ipynb` - Reproduce figure 3a in the manuscript.

### Tables
`./tables/Interpatient_SNVs_v2.csv.gz` - All SNV data from clinical samples (Aug 23rd, 2020).

`./tables/Wastewater_SNVs_v2.csv ` - All SNVs observed in wastewater samples

### Raw data
`./data/wastewater/` - Filtered bams and inStrain outputs of all wastewater samples.

`./data/viral_abundance_analysis/` - All data related to mapping to all RefSeq viruses for relative abundances.


## 2. Basic analysis pipeline

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
