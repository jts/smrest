def get_dir(wildcards):
    d = { "PAU59949.d052sup4305mCG_5hmCGvHg38": "colo829",
          "PAU61426.d052sup4305mCG_5hmCGvHg38": "colo829",
          "PAU59807.d052sup4305mCG_5hmCGvHg38": "colo829bl",
          "PAU61427.d052sup4305mCG_5hmCGvHg38": "colo829bl" }
    return d[wildcards.flowcell]

def get_sample_rate(wildcards):
    d = { "PAU59949.d052sup4305mCG_5hmCGvHg38": "0.277",
          "PAU61426.d052sup4305mCG_5hmCGvHg38": "0.274",
          "PAU59807.d052sup4305mCG_5hmCGvHg38": "0.303",
          "PAU61427.d052sup4305mCG_5hmCGvHg38": "0.388" }
    return d[wildcards.flowcell]

rule get_flowcell_bam:
    output:
        temp("data/flowcells/{flowcell}.chr20.20x.bam")
    params:
        dir=get_dir,
        sample_rate=get_sample_rate
    shell:
        "samtools view -b -o {output} -s {params.sample_rate} s3://ont-open-data/colo829_2024.03/basecalls/{params.dir}/sup/{wildcards.flowcell}.bam chr20" 

rule get_bam:
    input:
        fc=expand("data/flowcells/{fc}.chr20.20x.bam", fc=[ "PAU59949.d052sup4305mCG_5hmCGvHg38", "PAU61426.d052sup4305mCG_5hmCGvHg38", "PAU59807.d052sup4305mCG_5hmCGvHg38", "PAU61427.d052sup4305mCG_5hmCGvHg38" ])
    output:
        "data/COLO829.mixture.chr20.bam"
    shell:
        "samtools merge -o {output} {input}"

rule prepare_demo:
    input:
        "data/COLO829.mixture.chr20.bam", "data/COLO829.mixture.chr20.bam.bai"
