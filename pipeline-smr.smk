import csv

configfile: "/.mounts/labs/simpsonlab/users/jsimpson/projects/nanopore/2022.06.27.smrest_dev/etc/snakemake_config.yaml"

genotyping_window_size = 10
phasing_window_size = 100
calling_window_size = 10

def get_calling_windows(win_size_mb):
    MB=1000000
    window_size = win_size_mb * MB
    fai = config["reference"] + ".fai"
    chromosomes = [ "chr" + str(i) for i in range(0, 23) ]
    chromosomes.append("chrX")
    chromosomes.append("chrY")
    windows = list()
    with open(fai) as f:
        for line in f:
            fields = line.rstrip().split()
            chrom_name = fields[0]
            chrom_len = int(fields[1])

            if chrom_name not in chromosomes:
                continue

            for i in range(1, chrom_len, window_size):
                end = min(chrom_len, i + window_size)
                windows.append( f"{chrom_name}:{i}-{end}" )
    return windows

#
# top level rule, based on sample names defined in config/command line
#
#rule all:
#    input:
#        expand("smrest_calls/{sample}/{sample}.whatshap_phased.called_regions.bed", sample=config['samples']),
#        expand("smrest_calls/{sample}/{sample}.whatshap_phased.calls.vcf", sample=config['samples'])


#
# hmftools suite for purity/ploidy/copy number inference
#
rule run_amber:
    input:
        bam="data/{sample}.bam",
        bai="data/{sample}.bam.bai"
    output:
        directory("hmftools/{sample}/amber")
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    threads: 8
    shell:
        """java -jar -Xmx8G {config[hmftools_amber_jar]} \
                -tumor {wildcards.sample} \
                -tumor_bam {input.bam} \
                -output_dir {output} \
                -validation_stringency LENIENT \
                -threads {threads} \
                -ref_genome_version V38 \
                -loci {config[hmftools_resources]}/HMFtools-Resources_dna_pipeline_v5_31_38_copy_number_GermlineHetPon.38.vcf"""
    
rule run_cobalt:
    input:
        bam="data/{sample}.bam",
        bai="data/{sample}.bam.bai"
    output:
        directory("hmftools/{sample}/cobalt")
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    threads: 8
    shell:
        """java -jar -Xmx8G {config[hmftools_cobalt_jar]} \
                -tumor {wildcards.sample} \
                -tumor_bam {input.bam} \
                -output_dir {output} \
                -validation_stringency LENIENT \
                -threads {threads} \
                -tumor_only_diploid_bed {config[hmftools_resources]}/HMFtools-Resources_dna_pipeline_v5_31_38_copy_number_DiploidRegions.38.bed\
                -gc_profile {config[hmftools_resources]}/HMFtools-Resources_dna_pipeline_v5_31_38_copy_number_GC_profile.1000bp.38.cnp"""


rule run_purple:
    input:
        amber="hmftools/{sample}/amber",
        cobalt="hmftools/{sample}/cobalt"
    output:
    	base=directory("hmftools/{sample}/purple/"),
    	purity="hmftools/{sample}/purple/{sample}.purple.purity.tsv",
    	circos="hmftools/{sample}/purple/circos/{sample}.circos.conf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    threads: 8
    shell:
        """java -jar -Xmx8G {config[hmftools_purple_jar]} \
                -tumor {wildcards.sample} \
                -output_dir {output.base} \
                -amber {input.amber} \
                -cobalt {input.cobalt} \
                -gc_profile {config[hmftools_resources]}/HMFtools-Resources_dna_pipeline_v5_31_38_copy_number_GC_profile.1000bp.38.cnp \
                -ref_genome {config[reference]} \
                -threads 8 \
                -ref_genome_version V38 \
                -ensembl_data_dir {config[hmftools_resources]}/ensembl_data"""

#
# Genotype gnomad SNPs
#
rule get_on_target_gnomad:
    output:
        "resources/genotype_sites.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    threads: 1
    shell:
        "bcftools filter -R {config[giab_hc_bed]} {config[gnomad_common]} > {output}"

rule genotype_regions:
    input:
        bam="data/{sample}.bam",
        bai="data/{sample}.bam.bai",
        vcf="resources/genotype_sites.vcf"
    params:
        memory_per_thread="32G",
        extra_cluster_opt=""
    output:
        "smrest_genotype/{sample}/per_region/{sample}.region.{region}.gnomad_genotype.vcf"
    shell:
        "{config[smrest]} genotype-hets -c {input.vcf} -r {wildcards.region} -g {config[reference]} {input.bam} > {output}"

rule merge_region_genotype:
    input:
        vcfs=expand("smrest_genotype/{{sample}}/per_region/{{sample}}.region.{r}.gnomad_genotype.vcf", r=get_calling_windows(genotyping_window_size))
    output:
        "smrest_genotype/{sample}/{sample}.gnomad_genotype.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "bcftools concat {input.vcfs} | bcftools sort > {output}"

#
# Phase the genotyping results
#

rule compress_vcf:
    input:
        vcf="{prefix}.vcf"
    output:
        vcf="{prefix}.vcf.gz"
    threads: 1
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "bgzip {input}"

rule tabix_vcf:
    input:
        vcf="{prefix}.vcf.gz"
    output:
        tbi="{prefix}.vcf.gz.tbi"
    threads: 1
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "tabix {input}"

rule get_region_vcf:
    input:
        vcf="smrest_genotype/{sample}/{sample}.gnomad_genotype.vcf.gz",
        tbi="smrest_genotype/{sample}/{sample}.gnomad_genotype.vcf.gz.tbi"
    output:       
        vcf="phasing/{sample}/region_vcf/{sample}.gnomad_genotype.{region}.vcf"
    threads: 1
    params:
        memory_per_thread="8G",
        extra_cluster_opt=""
    shell:
        "bcftools view {input.vcf} {wildcards.region} > {output.vcf}"

rule phase_variants_whatshap:
    input:
        vcf="phasing/{sample}/region_vcf/{sample}.gnomad_genotype.{region}.vcf",
        bam="data/{sample}.bam"
    output:
        vcf="phasing/{sample}/phased_region_vcf/{sample}.gnomad_genotype.{region}.whatshap_phased.vcf"
    threads: 1
    params:
        memory_per_thread="24G",
        extra_cluster_opt=""
    shell:   
        "whatshap phase --ignore-read-groups -r {config[reference]} -o {output.vcf} {input.vcf} {input.bam}"

rule merge_phased_vcfs:
    input:
        vcfs=expand("phasing/{{sample}}/phased_region_vcf/{{sample}}.gnomad_genotype.{r}.whatshap_phased.vcf", r=get_calling_windows(phasing_window_size))
    output:
        "phasing/{sample}/{sample}.gnomad_genotype.whatshap_phased.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "bcftools concat {input.vcfs} | bcftools sort > {output}"

rule haplotag_bam:
    input:
        vcf="phasing/{sample}/{sample}.gnomad_genotype.whatshap_phased.vcf",
        bam="data/{sample}.bam"
    output:
        "haplotag/{sample}/{sample}.gnomad_genotype.whatshap_phased.tagged.bam"
    threads: 8
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "../longphase/longphase haplotag --log -m 1 -p 0.85 -t 8 -s {input.vcf} -b {input.bam} -o haplotag/{wildcards.sample}/{wildcards.sample}.gnomad_genotype.whatshap_phased.tagged"

def get_fai(wildcards):
    return config['reference'] + ".fai"

#
# Somatic mutation calling
#
def get_purity_from_purple(wildcards):
    fn = f"hmftools/{wildcards.sample}/purple/{wildcards.sample}.purple.purity.tsv"
    with open(fn) as f:
        rd = csv.DictReader(f, delimiter="\t", quotechar='"')
        for row in rd:
            p = float(row['purity'])
            return p

rule call_regions:
    input:
        bam="haplotag/{sample}/{sample}.gnomad_genotype.whatshap_phased.tagged.bam",
        bai="haplotag/{sample}/{sample}.gnomad_genotype.whatshap_phased.tagged.bam.bai"
    params:
        memory_per_thread="24G",
        extra_cluster_opt="",
        purity=get_purity_from_purple
    output:
        vcf="smrest_calls_phased_bam/{sample}/per_region/{sample}.whatshap_phased.region.{region}.vcf",
        bed="smrest_calls_phased_bam/{sample}/per_region/{sample}.whatshap_phased.region.{region}.bed"
    shell:
        "{config[smrest]} call -m haplotype-likelihood --purity {params.purity} -r {wildcards.region} -g {config[reference]} -o {output.bed} {input.bam} > {output.vcf}"

rule call_regions_internal_tag:
    input:
        #bam="haplotag/{sample}/{sample}.gnomad_genotype.whatshap_phased.tagged.bam",
        #bai="haplotag/{sample}/{sample}.gnomad_genotype.whatshap_phased.tagged.bam.bai"
        bam="data/{sample}.bam",
        bai="data/{sample}.bam.bai",
        phased_vcf="phasing/{sample}/{sample}.gnomad_genotype.whatshap_phased.vcf",
    	purity="hmftools/{sample}/purple/{sample}.purple.purity.tsv"
    params:
        memory_per_thread="24G",
        extra_cluster_opt="",
        purity=get_purity_from_purple
    output:
        vcf="smrest_calls/{sample}/per_region/{sample}.whatshap_phased.region.{region}.vcf",
        bed="smrest_calls/{sample}/per_region/{sample}.whatshap_phased.region.{region}.bed"
    shell:
        "{config[smrest]} call -m haplotype-likelihood --purity 0.5 -r {wildcards.region} -g {config[reference]} -p {input.phased_vcf} -o {output.bed} {input.bam} > {output.vcf}"

rule merge_region_calls:
    input:
        vcfs=expand("smrest_calls/{{sample}}/per_region/{{sample}}.whatshap_phased.region.{r}.vcf", r=get_calling_windows(calling_window_size))
    output:
        "smrest_calls/{sample}/{sample}.{phaser}.calls.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "bcftools concat {input.vcfs} | bcftools sort > {output}"

rule merge_region_bed:
    input:
        beds=expand("smrest_calls/{{sample}}/per_region/{{sample}}.whatshap_phased.region.{r}.bed", r=get_calling_windows(calling_window_size))
    output:
        "smrest_calls/{sample}/{sample}.{phaser}.called_regions.bed"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "cat {input.beds} > {output}"

rule final_bed:
    input:
        "smrest_calls/{sample}/{sample}.{phaser}.called_regions.bed"
    output:
        "smrest_calls/{sample}/{sample}.{phaser}.final_call_regions.bed"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """bedtools intersect -b {config[giab_hc_bed]} \
                              -a {input} > {output}"""

rule final_calls:
    input:
        vcf="smrest_calls/{sample}/{sample}.{phaser}.calls.vcf.gz",
        tbi="smrest_calls/{sample}/{sample}.{phaser}.calls.vcf.gz.tbi",
        bed="smrest_calls/{sample}/{sample}.{phaser}.final_call_regions.bed"
    output:
        "smrest_calls/{sample}/{sample}.{phaser}.final_q{min_somatic_qual}_pass_calls.vcf"
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        """
        if [ -s {input.bed} ]; then
            bcftools filter -T {input.bed} -i 'QUAL > {wildcards.min_somatic_qual} && (FILTER == \"PASS\" || FILTER == \".\")' {input.vcf} > {output}
        else
            zgrep \"^#\" {input.vcf} > {output}
        fi
        """

rule index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    threads: 1
    params:
        memory_per_thread="4G",
        extra_cluster_opt=""
    shell:
        "samtools index {input}"
