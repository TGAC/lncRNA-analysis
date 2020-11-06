# lncRNA-analysis
This is a Snakemake-based pipeline to annotate lncRNA using existing annotation. 

Briefly, it uses HISAT2 to align all reads to the transcriptome, then it discovers new genes using Stringtie, it creates new gene list with all known and novel genes. Then, it separates the novel genes and uses CPAT and CPC to annotate novel genes with coding potential and extract non-coding genes. It also identifies and removes transcripts where an isoform is coding.
