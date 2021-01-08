## Snakemake - lncRNA analysis
##
## Pipeline is to find non-coding RNAs from raw reads, 
## using reference genome and refenrece gene models.
##
## Anil S. Thanki
## GitHub: https://github.com/anilthanki/
## Anil.Thanki@earlham.ac.uk
##
##


## Configuration file containing absolute file locations
configfile: "ncRNA.json"

SAMPLES = ["ARI_D", "ARI_F", "ARI_G", "ARI_R", "ARI_S", "ARI_V"]

READS = ["1", "2"]


# --- Help Rules --- #
## help :                                       Prints help comments for Snakefile
rule help:
    input: "Snakefile"
    shell:
        "sed -n 's/^##//p' {input}"
        

## all :                                        List out final output files
rule all:
		input : expand([config["dir"]["outdir"]+'/FastQC/{sample}_{read}/{sample}_{read}_val_{read}_fastqc.html',  config["dir"]["outdir"]+"/GTF_nc_transcript_filtered/{sample}.gtf"], sample=SAMPLES, read=READS) 
		


## trimgalore :                                 Trims reads and removes adapters
rule trimgalore:
	input:
		fq1 = config["dir"]["raw_datadir"]+"/{sample}_1.fq.gz",
		fq2 = config["dir"]["raw_datadir"]+"/{sample}_2.fq.gz",
	output:
		fq1 = temp(config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_1_val_1.fq.gz"),  
		fq2 = temp(config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_2_val_2.fq.gz"),
	params:
		dir = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}",
	conda:
	    "conda/trimgalore.yaml"
	shell:
		"""
			trim_galore --illumina -o {params.dir} --paired {input.fq1} {input.fq2}
		"""


## qc :                                         Performs quality check on trimmed reads
rule qc:
	input:
		fq = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_{read}_val_{read}.fq.gz",
	output:
		fq = config["dir"]["outdir"]+"/FastQC/{sample}_{read}/{sample}_{read}_val_{read}_fastqc.html",
	params:
		fq_dir = config["dir"]["outdir"]+"/FastQC/{sample}_{read}/"
	conda:
	    "conda/fastqc.yaml"
	shell:
		"""
			mkdir -p {params.fq_dir} && \
			fastqc {input.fq} --outdir={params.fq_dir}
		"""


## align_hista :                                Performs alignments using HISAT on trimmed reads
##                                              SAM files are temp files, which will be deleted after all rules that use it as an input are completed
rule align_hisat:
	input:
		fq1 = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_1_val_1.fq.gz",
		fq2 = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_2_val_2.fq.gz",
	output: 
		sam = temp(config["dir"]["outdir"]+"/SAM/{sample}.sam"),
		summary = config["dir"]["outdir"]+"/SAM/{sample}_summary.txt",
	threads: 4
	params:
	    index = lambda wc: config['hisat_index'][wc.sample],
	wildcard_constraints:
	    sample = '\w+'
	conda:
	    "conda/hisat.yaml"
	shell:
		"""
			hisat2 -q -p {threads} --dta -x {params.index} -1 {input.fq1} -2 {input.fq2} -S {output.sam} --summary {output.summary} 
		"""


## sam2bam :                                    Converts SAM to BAM format
##                                              BAM files are temp files, which will be deleted after all rules that use it as an input are completed
rule sam2bam:
	input:  config["dir"]["outdir"]+"/SAM/{sample}.sam"
	output: temp(config["dir"]["outdir"]+"/BAM/{sample}.bam")
	conda:
	    "conda/samtools.yaml"
	shell:
		"""
			samtools view -b {input} > {output} 
		"""


## sort_bam :                                   Sorts BAM files
rule sort_bam:
	input:  config["dir"]["outdir"]+"/BAM/{sample}.bam"
	output: config["dir"]["outdir"]+"/BAM/{sample}_sorted.bam"
	conda:
	    "conda/samtools.yaml"
	shell:
		"""
			samtools sort {input} > {output}
		"""


## index_bam :                                  Indexes BAM alignments
rule index_bam:
	input:	config["dir"]["outdir"]+"/BAM/{sample}_sorted.bam"
	output: config["dir"]["outdir"]+"/BAM/{sample}.bam.csi"
	conda:
	    "conda/samtools.yaml"
	shell:
		"""
			samtools index -c {input} {output}
		"""


## stringtie :                                  StringTie is used to annotate gene models from the alignment files using provided gene information as reference
##                                              Here also produces abundance files
rule stringtie:
	input: 
		bam = config["dir"]["outdir"]+"/BAM/{sample}_sorted.bam",
		gff = lambda wc: config['sample_gff'][wc.sample]
	output:
		gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_stringtie.gtf",
		abundance = config["dir"]["outdir"]+"/GTF_stringtie/{sample}.out"
	params:
		ballgown = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_ballgown_input/",
		label = "{sample}-Stringtie"
	wildcard_constraints:
	    sample = '\w+'
	conda:
	    "conda/stringtie.yaml"
	shell:
		"""
			stringtie {input.bam} -p 4 -G {input.gff} -o {output.gtf} -A {output.abundance} -b {params.ballgown} -l {params.label}
		"""


## gffcompare_for_novel_transcripts :           GFFcompare is used to find novel transcript annotated by StringTie
##                                              Class 'u' is used to fetch novel genes
rule gffcompare_for_novel_transcripts:
	input:
		ref_gff = lambda wc: config['sample_gff'][wc.sample],
		query_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_stringtie.gtf",
	output:
		novel_transcripts_list = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.txt",
		novel_transcripts_list2 = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts_with_quotes.txt",
		novel_transcripts = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.gtf",
	params:
		gtf_dir = config["dir"]["outdir"]+"/GTF_stringtie/",
		label = "{sample}",
		gffcompare_file = "{sample}.{sample}_stringtie.gtf.tmap",
		column = "-f5" # -f4 for genes
	wildcard_constraints:
	    sample = '\w+'
	conda:
	    "conda/gffcompare.yaml"
	shell:
		"""
			cd {params.gtf_dir} && \
			gffcompare -r {input.ref_gff} {input.query_gtf} -o {params.label} && \
			cat {params.gffcompare_file} | awk '$3=="u"{{print $0}}' | cut {params.column} | sort | uniq > {output.novel_transcripts_list} && \
			sed 's/^/"/; s/$/"/' {output.novel_transcripts_list} > {output.novel_transcripts_list2} && \
			grep -F -f {output.novel_transcripts_list2} {input.query_gtf} > {output.novel_transcripts}	
		"""


## novel_transcript_fasta :                     GFFread is used to fetch FASTA for novel transcript from reference genome
rule novel_transcripts_fasta:
	input:
		novel_transcripts_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.gtf",
		reference_fasta = lambda wc: config['ref_genome'][wc.sample],
		#reference_fasta = config["dir"]["ref_genome"]
	output:
		novel_transcripts_fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa"
	wildcard_constraints:
	    sample = '\w+'
	conda:
	    "conda/gffread.yaml"
	shell:
		"""
			gffread -w {output.novel_transcripts_fasta} -g {input.reference_fasta} {input.novel_transcripts_gtf}
		"""


## cpat_training_data :                         Prepare hexamer files used in CPAT run using provided training dataset
rule CPAT_training_data:
	input:
		cds = config["dir"]["cpat"]["training-data"]["cds"],
		ncrna = config["dir"]["cpat"]["training-data"]["ncrna"]
	output:
		hexamer = config["dir"]["outdir"]+"/cpat/models/wheat.tsv"
	conda:
	    "conda/cpat.yaml"
	shell:
		"""
			make_hexamer_tab.py -c {input.cds} -n {input.ncrna} > {output.hexamer}
		"""


## cpat_trainig_model :                         Prepare logitmodel files used in CPAT run using provided training dataset
rule cpat_trainig_model:
	input:
		cds = config["dir"]["cpat"]["training-data"]["cds"],
		ncrna = config["dir"]["cpat"]["training-data"]["ncrna"],
		hexamer = config["dir"]["outdir"]+"/cpat/models/wheat.tsv"
	params:
		logitmodel = config["dir"]["outdir"]+"/cpat/models/wheat",
	output:
		RData = config["dir"]["outdir"]+"/cpat/models/wheat.logit.RData"
	conda:
	    "conda/cpat.yaml"
	shell:
		"""
			make_logitModel.py -c {input.cds} -n {input.ncrna} -x {input.hexamer} -o {params.logitmodel}
		"""


## cpat :                                       Predict coding potential of novel transcripts using CPAT
rule cpat:
	input:
		fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa",
		logitmodel = config["dir"]["outdir"]+"/cpat/models/wheat.logit.RData",
		hexamer = config["dir"]["outdir"]+"/cpat/models/wheat.tsv"
	output:
		tsv = config["dir"]["outdir"]+"/cpat/output/{sample}.tsv",
		tsv_lowercase = config["dir"]["outdir"]+"/cpat/output/{sample}.lowercase.tsv",
	conda:
	    "conda/cpat.yaml"
	shell:
		"""
			cpat.py -g {input.fasta}  -d {input.logitmodel} -x {input.hexamer} -o {output.tsv} && \
			
			# changing the case of Stringtie, CPAT changes Stringtie to STRINGTIE in previous rule 
			perl -p -e 's/STRINGTIE/Stringtie/g' {output.tsv} > {output.tsv_lowercase}
		"""
		


## filter_cpat :                                Separates non-coding transcripts from CPAT analysis using coding potential cutoff
rule filter_cpat:
	input:
		tsv = config["dir"]["outdir"]+"/cpat/output/{sample}.lowercase.tsv",
	output:
		nctsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.tsv",
	params:
		cutoff = 0.4
	shell:
		"""
			awk -F "\t" '{{ if($6 <= {params.cutoff} ) {{print $1}} }}' {input.tsv} > {output.nctsv}
		"""


## sort_cpat :                                  Sorts non-coding transcripts for comparison with CPC prediction
rule sort_cpat:
	input:
		nctsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.tsv",
	output:
	    ncsortedtsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.sorted.tsv",
	shell:
		"""
			sort {input.nctsv} > {output.ncsortedtsv}
		"""




## cpc :                                        Predict coding potential of novel transcripts using CPC
rule cpc:
	input:
		fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa",
	output:
		tsv = config["dir"]["outdir"]+"/cpc/output/{sample}.tsv",
	conda:
	    "conda/cpc.yaml"
	shell:
		"""
			CPC2.py -i {input.fasta}  -o {output.tsv}
		"""


## filter_cpc :                                 Separates non-coding transcripts from CPC analysis
rule filter_cpc:
	input:
		tsv = config["dir"]["outdir"]+"/cpc/output/{sample}.tsv",
	output:
		nctsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.tsv",
	params:
		filter = "noncoding"	
	shell:
		"""
			awk '{{ if($8 == "{params.filter}" ) {{print $1}} }}' {input.tsv} > {output.nctsv}
		"""


## sort_cpc :                                   Sorts non-coding transcripts for comparison with CPAT prediction
rule sort_cpc:
	input:
		nctsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.tsv",
	output:
		ncsortedtsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.sorted.tsv",
	shell:
		"""
			sort {input.nctsv} > {output.ncsortedtsv}
		"""


## cpat_cpc_intersect :                         Finds common non-coding transcripts predicted by both CPAT and CPC
rule cpat_cpc_intersect:
    input:
        cpat = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.sorted.tsv",
        cpc = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.sorted.tsv",
    output:
        comm = config["dir"]["outdir"]+"/nc/{sample}.common_nc.tsv",
    shell:
        """
            comm -12 {input.cpat} {input.cpc} > {output.comm}
		"""


## non_coding_transcript :                      Fetches non-coding transcripts using output of 'cpat_cpc_intersect' from novel transcripts generated by StringTie
rule non_coding_transcript_from_cpat:
	input:
		tsv =  config["dir"]["outdir"]+"/nc/{sample}.common_nc.tsv",
		query_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.gtf",
	output:
		tsv =  config["dir"]["outdir"]+"/nc/{sample}.quoted_common_nc.tsv",
		novel_transcripts = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
	shell:
		"""
			sed 's/^/"/; s/$/"/' {input.tsv} > {output.tsv} && \
			grep -F -f {output.tsv} {input.query_gtf} > {output.novel_transcripts}
		"""


## genelist_for_non_coding_transcripts :        Fetches gene list from GTF of 'non_coding_transcript'
rule genelist_for_non_coding_transcripts:
	input:
		nc_gtf =  config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
	output:
		gene_list = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.tsv"
	shell:
		"""
			#keeping " to avoid matching with wrong gene ids
			perl -lne 'print @m if @m=(/((?:gene_id)\s+\S+)/g);' {input.nc_gtf} | perl -p -i -e 's/gene_id|;| //g' > {output.gene_list} && \
			sort {output.gene_list} | uniq
		"""



## all_genes_for_non_coding_transcript :        Fetches all transcript of gene list to look for coding isoforms
rule all_genes_for_non_coding_transcript:
	input:
		query_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_stringtie.gtf",
		gene_list = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.tsv"
	output:
		transcripts = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.gtf"
	shell:
		"""
			# gene list is already quoted so not adding quotes again, could use Mikado
			grep -F -f {input.gene_list} {input.query_gtf} > {output.transcripts}
		"""


## gene_transcript_pair :                       Prepares gene id, transcript id list from GTF of 'non_coding_transcript' and 'all_genes_for_non_coding_transcript'
##                                              First GTF contains non-coding transcripts and later contains all isoforms of non-codong transcripts
rule gene_transcript_pair:
	input:
		nc_gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
		combine_gtf = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.gtf",
	output:
		nc_pair = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.pair.tsv",
		combine_pair = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.pair.tsv",
	params: 
			scripts = config["scripts"]
	shell:
		"""
			python {params.scripts}/tra.py < {input.nc_gtf} > {output.nc_pair} && \
			python {params.scripts}/tra.py < {input.combine_gtf} > {output.combine_pair}
		"""


## coding_isoform_of_non_coding_transcripts :   Prepares a list of genes which has coding isoforms of non-codong transcripts
rule coding_isoform_of_non_coding_transcripts:
	input:
		nc_pair = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.pair.tsv",
		combine_pair = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.pair.tsv",
	output:
		diff = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.ncgenes.tsv",
	params: 
			scripts = config["scripts"]
	shell:
		"""
			python {params.scripts}/filter.py {input.nc_pair} {input.combine_pair} > {output.diff}
		"""


## non_coding_genes :                           Prepares GTF of the genes with only non-codong transcripts
##                                              Adds quotes around id using sed 
##                                              uses reverse grep to fetch genes which are not present in input list 
rule non_coding_genes:
    input:
        gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
        genes = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.ncgenes.tsv"
    output:
        quoted_genes = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.ncgenes.quoted.tsv",
        gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.filtered.gtf"
    shell:
        """
            sed 's/^/"/; s/$/"/' {input.genes} > {output.quoted_genes} && \
            grep -v -F -f {output.quoted_genes} {input.gtf} > {output.gtf}
        """



## filter_non_coding_genes :                    Filters non-coding genes based on transcript length 
rule filter_non_coding_genes:
	input: config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.filtered.gtf",
	output: config["dir"]["outdir"]+"/GTF_nc_transcript_filtered/{sample}.gtf",
	conda:
	    "conda/gffread.yaml"
	shell:
		"""
			gffread {input} -T -l 200 -o {output}
		"""
