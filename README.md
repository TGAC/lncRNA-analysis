# lncRNA-analysis
A Snakemake-based pipeline to annotate lncRNA using existing annotation. 

## Description
This pipeline uses HISAT2 to align all reads to the transcriptome, then it discovers new genes using Stringtie, it creates new gene list with all known and novel genes. Then, it separates the novel genes and uses CPAT and CPC to annotate novel genes with coding potential and extract non-coding genes. It also identifies and removes transcripts where an isoform is coding.

## Pipeline
A simple overview of pipeline is as below:

<img src="img/pipeline.svg" height=1000px title="pipeline-overview" />

## Authors

### Developers:

- Wilfried Haerty (Earlham Institute)
- Anil S. Thanki (Earlham Institute)

