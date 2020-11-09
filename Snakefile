configfile: "ncRNA.json"

SAMPLES = ["ARI_D", "ARI_F", "ARI_G", "ARI_R", "ARI_S", "ARI_V"]
#SAMPLES = ["JAG_D", "JAG_F", "JAG_R", "JAG_S", "JAG_V"]
#SAMPLES = ["JUL_D", "JUL_F", "JUL_G", "JUL_R", "JUL_S", "JUL_V"]
#SAMPLES = ["LAN_D", "LAN_F", "LAN_G", "LAN_R", "LAN_S", "LAN_V"]
#SAMPLES = ["LER_D", "LER_F", "LER_R", "LER_S", "LER_V"]
#SAMPLES = ["MAC_D", "MAC_F", "MAC_R", "MAC_S", "MAC_V"]
#SAMPLES = ["MAT_D", "MAT_F", "MAT_G", "MAT_R", "MAT_S", "MAT_V"]
#SAMPLES = ["NOR_D", "NOR_F", "NOR_G", "NOR_R", "NOR_S", "NOR_V"]
#SAMPLES = ["STA_D", "STA_F", "STA_R", "STA_S", "STA_V"]

READS = ["1", "2"]
#, "ARI_G", "ARI_R", "ARI_S", "ARI_V"]


# Define input files
def cpat_files(wildcards):
	files = expand( config["dir"]["outdir"]+"/cpat/output/{sample}_nc.tsv", sample=SAMPLES)
	return files

def get_gff(sample):
	return config['sample_gff'][sample]


rule all:
		input : expand([config["dir"]["outdir"]+'/FastQC/{sample}_{read}/{sample}_{read}_val_{read}_fastqc.html',  config["dir"]["outdir"]+"/GTF_nc_transcript_filtered/{sample}.gtf"], sample=SAMPLES, read=READS) 


rule trimgalore:
	input:
		fq1 = config["dir"]["raw_datadir"]+"/{sample}_1.fq.gz",
		fq2 = config["dir"]["raw_datadir"]+"/{sample}_2.fq.gz",
	output:
		fq1 = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_1_val_1.fq.gz",  
				fq2 = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_2_val_2.fq.gz",
	params:
		dir = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}",
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/trimgalore && \
			trim_galore --illumina -o {params.dir} --paired {input.fq1} {input.fq2}
		"""

rule qc:
	input:
		fq = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_{read}_val_{read}.fq.gz",
	output:
		fq = config["dir"]["outdir"]+"/FastQC/{sample}_{read}/{sample}_{read}_val_{read}_fastqc.html",
	params:
		fq_dir = config["dir"]["outdir"]+"/FastQC/{sample}_{read}/"	
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/fastqc && \
			mkdir -p {params.fq_dir} && \
			fastqc {input.fq} --outdir={params.fq_dir}
		"""

rule align_hisat:
	input:
		fq1 = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_1_val_1.fq.gz",
		fq2 = config["dir"]["outdir"]+"/Trimmed_Reads/{sample}/{sample}_2_val_2.fq.gz",
	output: 
		sam = temp(config["dir"]["outdir"]+"/SAM/{sample}.sam"),
		summary = config["dir"]["outdir"]+"/SAM/{sample}_summary.txt",
	threads: 4
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/hisat && \
			hisat2 -q -p {threads} --dta -x /ei/projects/9/90a2cc84-e6dd-40a3-ab23-38dcd2389ff0/data/reference/Triticum_aestivum.IWGSC.dna.toplevel -1 {input.fq1} -2 {input.fq2} -S {output.sam} --summary {output.summary} 
		"""

rule sam2bam:
	input:  config["dir"]["outdir"]+"/SAM/{sample}.sam"
	output: temp(config["dir"]["outdir"]+"/BAM/{sample}.bam")
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/samtools && \
			samtools view -b {input} > {output} 
		"""

rule sort_bam:
	input:  config["dir"]["outdir"]+"/BAM/{sample}.bam"
	output: config["dir"]["outdir"]+"/BAM/{sample}_sorted.bam"
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/samtools && \
			samtools sort {input} > {output}
		"""

rule index_bam:
	input:	config["dir"]["outdir"]+"/BAM/{sample}_sorted.bam"
	output: config["dir"]["outdir"]+"/BAM/{sample}.bam.csi"
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/samtools && \
			samtools index -c {input} {output}
		"""


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
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/stringtie && \
			stringtie {input.bam} -p 4 -G {input.gff} -o {output.gtf} -A {output.abundance} -b {params.ballgown} -l {params.label}
		"""

rule gffcompare_for_novel_transcripts:
	input:
#		ref_gff = lambda wc: get_gff('{sample}'.format(sample=wc.sample)),
		ref_gff = lambda wc: config['sample_gff'][wc.sample],
		query_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_stringtie.gtf",
	output:
		novel_transcripts_list = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.txt",
		novel_transcripts = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.gtf",
	params:
		gtf_dir = config["dir"]["outdir"]+"/GTF_stringtie/",
		label = "{sample}",
		gffcompare_file = "{sample}.{sample}_stringtie.gtf.tmap",
		column = "-f5" # -f4 for genes
	wildcard_constraints:
	    sample = '\w+'
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/gffcompare && \
			cd {params.gtf_dir} && \
			gffcompare -r {input.ref_gff} {input.query_gtf} -o {params.label} && \
			cat {params.gffcompare_file} | awk '$3=="u"{{print $0}}' | cut {params.column} | sort | uniq > {output.novel_transcripts_list} && \
			grep -F -f {output.novel_transcripts_list} {input.query_gtf} > {output.novel_transcripts}	
		"""

rule novel_transcripts_fasta:
	input:
		novel_transcripts_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.gtf",
		reference_fasta = config["dir"]["ref_genome"]
	output:
		novel_transcripts_fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa"
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/gffread && \
			gffread -w {output.novel_transcripts_fasta} -g {input.reference_fasta} {input.novel_transcripts_gtf}
		"""

rule CPAT_training_data:
	input:
		cds = config["dir"]["cpat"]["training-data"]["cds"],
		ncrna = config["dir"]["cpat"]["training-data"]["ncrna"]
	output:
		hexamer = config["dir"]["outdir"]+"/cpat/models/wheat.tsv"
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/cpat && \
			make_hexamer_tab.py -c {input.cds} -n {input.ncrna} > {output.hexamer}
		"""

rule cpat_trainig_model:
	input:
		cds = config["dir"]["cpat"]["training-data"]["cds"],
		ncrna = config["dir"]["cpat"]["training-data"]["ncrna"],
		hexamer = config["dir"]["outdir"]+"/cpat/models/wheat.tsv"
	params:
		logitmodel = config["dir"]["outdir"]+"/cpat/models/wheat",
	output:
		RData = config["dir"]["outdir"]+"/cpat/models/wheat.logit.RData"
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/cpat && \
			make_logitModel.py -c {input.cds} -n {input.ncrna} -x {input.hexamer} -o {params.logitmodel}
		"""

rule cpat:
	input:
		fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa",
		logitmodel = config["dir"]["outdir"]+"/cpat/models/wheat.logit.RData",
		hexamer = config["dir"]["outdir"]+"/cpat/models/wheat.tsv"
	output:
		tsv = config["dir"]["outdir"]+"/cpat/output/{sample}.tsv",
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/cpat && \
			cpat.py -g {input.fasta}  -d {input.logitmodel} -x {input.hexamer} -o {output.tsv} && \
			
			# changing the case of Stringtie, CPAT changes Stringtie to STRINGTIE in previous rule 
			perl -p -i -e 's/STRINGTIE/Stringtie/g' {output.tsv}
		"""
		
rule filter_cpat:
	input:
		tsv = config["dir"]["outdir"]+"/cpat/output/{sample}.tsv",
	output:
		nctsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.tsv",
	params:
		cutoff = 0.4
	shell:
		"""
			awk -F "\t" '{{ if($6 <= {params.cutoff} ) {{print $1}} }}' {input.tsv} > {output.nctsv}
		"""
		
rule sort_cpat:
	input:
		nctsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.tsv",
	output:
	    ncsortedtsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc_sorted.tsv",
	shell:
		"""
			sort {input.nctsv} > {output.ncsortedtsv}
		"""

rule cpc:
	input:
		fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa",
	output:
		tsv = config["dir"]["outdir"]+"/cpc/output/{sample}.tsv",
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/cpc && \
			CPC2.py -i {input.fasta}  -o {output.tsv}
		"""
		
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

rule sort_cpc:
	input:
		nctsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.tsv",
	output:
		ncsortedtsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc_sorted.tsv",
	shell:
		"""
			sort {input.nctsv} > {output.ncsortedtsv}
		"""
		
rule cpat_cpc_intersect:
    input:
        cpat = config["dir"]["outdir"]+"/cpat/output/{sample}.nc_sorted.tsv",
        cpc = config["dir"]["outdir"]+"/cpc/output/{sample}.nc_sorted.tsv",
    output:
        comm = config["dir"]["outdir"]+"/nc/{sample}.common_nc.tsv",
    shell:
        """
            comm -12 {input.cpat} {input.cpc} > {output.comm}
		"""
		
rule non_coding_transcript_from_cpat:
	input:
		#for cpc
		tsv =  config["dir"]["outdir"]+"/nc/{sample}.common_nc.tsv",
		#tsv =  config["dir"]["outdir"]+"/cpat/output/{sample}_nc_sorted.tsv",
		query_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_novel_transcripts.gtf",
	output:
		novel_transcripts = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
	shell:
		"""
			grep -F -f {input.tsv} {input.query_gtf} > {output.novel_transcripts}
		"""

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
			
rule all_genes_for_non_coding_transcript:
	input:
		query_gtf = config["dir"]["outdir"]+"/GTF_stringtie/{sample}_stringtie.gtf",
		gene_list = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.tsv"
	output:
		transcripts = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.gtf"
	shell:
		"""
			grep -F -f {input.gene_list} {input.query_gtf} > {output.transcripts}
		"""
			
rule gene_transcript_pair:
	input:
		nc_gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
		combine_gtf = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.gtf",
	output:
		nc_pair = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.pair.tsv",
		combine_pair = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.pair.tsv",
	shell:
		"""
			python /hpc-home/thankia/MyComputer/scripts/test/tra.py < {input.nc_gtf} > {output.nc_pair} && \
			python /hpc-home/thankia/MyComputer/scripts/test/tra.py < {input.combine_gtf} > {output.combine_pair}
		"""

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

rule non_coding_genes:
    input:
        gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
        genes = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.ncgenes.tsv"
    output:
        gff = temp(config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gff3"),
        gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.filtered.gtf",
        filtered_gff = temp(config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.filtered.gff3")
    shell:
        """
            source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/mikado && \
            mikado util convert -if gtf -of gff3 {input.gtf} > {output.gff} && \
            mikado util grep -v --genes {input.genes} {output.gff} > {output.filtered_gff} && \
            mikado util convert -if gff3 -of gtf {output.filtered_gff} > {output.gtf}
        """


rule filter_non_coding_genes:
	input: config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.filtered.gtf",
	output: config["dir"]["outdir"]+"/GTF_nc_transcript_filtered/{sample}.gtf",
	shell:
		"""
			source activate /ei/software/testing/python_miniconda/4.5.4_py3.6_cs/x86_64/envs/gffcompare && \
			gffread {input} -T -l 200 -o {output}
		"""
