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

