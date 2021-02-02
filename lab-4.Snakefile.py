# default rule that tells snakemake to create the .vcf file
# if it is not run with any specific rule or file request.

rule all:
    input:
        "SRR2584857_1.ecoli-rel606.vcf",

rule download_data:
    conda: "env-wget.yml"
    output: "SRR2584857_1.fastq.gz"
    shell: """
        wget https://osf.io/4rdza/download -O {output}
    """

rule download_genome:
    conda: "env-wget.yml"
    output: "ecoli-rel606.fa.gz"
    shell:
        "wget https://osf.io/8sm92/download -O {output}"

rule map_reads:
    conda: "env-minimap.yml"
    input:
        ref="ecoli-rel606.fa.gz",
        reads="{sample}.fastq.gz"
    output: "{sample}.ecoli-rel606.sam"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """

rule sam_to_bam:
    message: "Convert text SAM file to binary BAM file for sample {wildcards.sample}"
    conda: "env-minimap.yml"
    input: "{sample}.ecoli-rel606.sam"
    output: "{sample}.ecoli-rel606.bam"
    shell: """
        samtools view -b -F 4 {input} > {output}
     """

rule sort_bam:
    conda: "env-minimap.yml"
    input: "{sample}.ecoli-rel606.bam"
    output: "{sample}.ecoli-rel606.bam.sorted"
    shell: """
        samtools sort {input} > {output}
    """

rule call_variants:
    conda: "env-bcftools.yml"
    input:
        ref="ecoli-rel606.fa.gz",
        bamsort="{sample}.ecoli-rel606.bam.sorted"
    output:
        refout = "{sample}.ecoli-rel606.fa",
        pileup="{sample}.ecoli-rel606.pileup",
        bcf="{sample}.ecoli-rel606.bcf",
        vcf="{sample}.ecoli-rel606.vcf"
    shell: """
        gunzip -c {input.ref} > {output.refout}
        bcftools mpileup -Ou -f {output.refout} {input.bamsort} > {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view -i '%QUAL>=20' {output.bcf} > {output.vcf}
    """
