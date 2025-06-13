#!/usr/bin/env python3

base = "data/"
suffix_1 = "_1.fastq"
suffix_2 = "_2.fastq"

# Define sample IDs using glob_wildcards
SAMPLES, = glob_wildcards(base + "{sample}" + suffix_1)

# Define all desired outputs
want_all = []
want_all.append(expand("data/{sample}_1.fastq", sample=SAMPLES))
want_all.append(expand("data/{sample}_2.fastq", sample=SAMPLES))
want_all.append("references/GCF_001617525.2_ASM161752v2_genomic.fna")
want_all.append("references/GCF_001617525.2_ASM161752v2_genomic.gff")
want_all.append("references/GCF_001617525.2.gtf")
want_all.append("references/genome_index.1.ht2")
want_all.append(expand("data/fastqc_results/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["1", "2"]))
want_all.append("data/fastqc_results/multiqc_report.html")
want_all.append(expand("results/{sample}.bam", sample=SAMPLES))
want_all.append("results/counts.txt")
want_all.append("results/simple_counts.txt")
want_all.append("Tabel_DGE.csv")

rule all:
    input: want_all

rule unzip_fastq:
    input:
        read = "data/{sample}_{read}.fastq.gz"
    output:
        read = "data/{sample}_{read}.fastq"
    threads: 1
    shell: """
        gunzip -c {input.read} > {output.read}
    """

rule download_reference:
    output:
        ref = "references/GCF_001617525.2_ASM161752v2_genomic.fna.gz"
    threads: 1
    shell: """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/617/525/GCF_001617525.2_ASM161752v2/GCF_001617525.2_ASM161752v2_genomic.fna.gz -O {output.ref}
    """

rule unzip_reference:
    input:
        ref = "references/GCF_001617525.2_ASM161752v2_genomic.fna.gz"
    output:
        ref = "references/GCF_001617525.2_ASM161752v2_genomic.fna"
    threads: 1
    shell: """
        gunzip -c {input.ref} > {output.ref}
    """

rule download_gff:
    output:
        gff = "references/GCF_001617525.2_ASM161752v2_genomic.gff.gz"
    threads: 1
    shell: """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/617/525/GCF_001617525.2_ASM161752v2/GCF_001617525.2_ASM161752v2_genomic.gff.gz -O {output.gff}
    """

rule unzip_gff:
    input:
        gff = "references/GCF_001617525.2_ASM161752v2_genomic.gff.gz"
    output:
        gff = "references/GCF_001617525.2_ASM161752v2_genomic.gff"
    threads: 1
    shell: """
        gunzip -c {input.gff} > {output.gff}
    """

rule hisat2_index:
    input:
        ref = "references/GCF_001617525.2_ASM161752v2_genomic.fna"
    output:
        idx = "references/genome_index.1.ht2"
    threads: 32
    shell: """
        hisat2-build {input.ref} references/genome_index
    """

rule fastqc:
    input:
        read = "data/{sample}_{read}.fastq"
    output:
        html = "data/fastqc_results/{sample}_{read}_fastqc.html",
        zip = "data/fastqc_results/{sample}_{read}_fastqc.zip"
    threads: 1
    shell: """
        fastqc {input.read} -o data/fastqc_results
    """

rule multiqc:
    input:
        expand("data/fastqc_results/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["1", "2"])
    output:
        report = "data/fastqc_results/multiqc_report.html"
    threads: 1
    shell: """
        multiqc data/fastqc_results/ -o data/fastqc_results
    """

rule hisat2_align:
    input:
        read1 = "data/{sample}_1.fastq",
        read2 = "data/{sample}_2.fastq",
        idx = "references/genome_index.1.ht2"
    output:
        bam = "results/{sample}.bam"
    threads: 32
    shell: """
        hisat2 -x references/genome_index -1 {input.read1} -2 {input.read2} | samtools sort > {output.bam}
    """

rule gff_to_gtf:
    input:
        gff = "references/GCF_001617525.2_ASM161752v2_genomic.gff"
    output:
        gtf = "references/GCF_001617525.2.gtf"
    threads: 1
    shell: """
        gffread {input.gff} -T -o {output.gtf}
    """

rule feature_counts:
    input:
        bam = expand("results/{sample}.bam", sample=SAMPLES),
        gtf = "references/GCF_001617525.2.gtf"
    output:
        counts = "results/counts.txt"
    threads: 32
    shell: """
        featureCounts -p --countReadPairs -t transcript -g transcript_id -a {input.gtf} -o {output.counts} {input.bam}
    """

rule simplify_counts:
    input:
        counts = "results/counts.txt"
    output:
        simple_counts = "results/simple_counts.txt"
    threads: 1
    shell: """
        cut -f 1,7-12 {input.counts} > {output.simple_counts}
    """

rule run_deseq2:
    input:
        counts = "results/simple_counts.txt",
        script = "deseq2.r"
    output:
        results = "Tabel_DGE.csv"
    threads: 1
    shell: """
        Rscript {input.script} {input.counts} {output.results}
    """
