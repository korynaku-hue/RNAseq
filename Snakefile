configfile: "config.yaml"

RAW_DIR = config["paths"]["raw_dir"]
OUT_DIR = config["paths"]["out_dir"]
THREADS = config["params"]["threads"]
SAMPLES = glob_wildcards(f"{RAW_DIR}/{{sample}}_1.fastq").sample
INDEX = multiext(f"{OUT_DIR}/04_index/genome",".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")

rule all:
    input:
        expand(f"{OUT_DIR}/01_qc_raw/{{sample}}_{{read}}_fastqc.html",sample=SAMPLES, read=[1, 2]),
        expand(f"{OUT_DIR}/03_qc_trimmed/{{sample}}.pair{{read}}.truncated_fastqc.html", sample=SAMPLES, read=[1, 2]),
        f"{OUT_DIR}/06_counts/gene_counts.txt"

#--READ PREPROCESSING--------------------------------------------------------

rule fastqc_raw:
    input:
        r1 = f"{RAW_DIR}/{{sample}}_1.fastq",
        r2 = f"{RAW_DIR}/{{sample}}_2.fastq"
    output:
        f"{OUT_DIR}/01_qc_raw/{{sample}}_1_fastqc.html",
        f"{OUT_DIR}/01_qc_raw/{{sample}}_2_fastqc.html"
    threads: THREADS
    conda: "env.yaml"
    shell:
        "fastqc -t {threads} -o {OUT_DIR}/01_qc_raw {input.r1} {input.r2}"

rule trim_reads:
    input:
        r1 = f"{RAW_DIR}/{{sample}}_1.fastq",
        r2 = f"{RAW_DIR}/{{sample}}_2.fastq"
    output:
        r1 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair1.truncated.gz",
        r2 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair2.truncated.gz"
    threads: THREADS
    conda: "env.yaml"
    shell:
        "AdapterRemoval --file1 {input.r1} --file2 {input.r2} "
        "--basename {OUT_DIR}/02_trimmed/{wildcards.sample} "
        "--threads {threads} --gzip"

rule fastqc_trimmed:
    input:
        r1 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair1.truncated.gz",
        r2 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair2.truncated.gz"
    output:
        f"{OUT_DIR}/03_qc_trimmed/{{sample}}.pair1.truncated_fastqc.html",
        f"{OUT_DIR}/03_qc_trimmed/{{sample}}.pair2.truncated_fastqc.html"
    threads: THREADS
    conda: "env.yaml"
    shell:
        "fastqc -t {threads} -o {OUT_DIR}/03_qc_trimmed {input.r1} {input.r2}"

#--ALIGNMENT-BASED QUANTIFICATION-------------------------------------------

rule prepare_ref:
    input:
        fa = config["paths"]["genome_gz"],
        gtf = config["paths"]["gtf_gz"]
    output:
        fa = f"{OUT_DIR}/00_reference/genome.fa",
        gtf = f"{OUT_DIR}/00_reference/annotation.gtf"
    shell:
        "gunzip -c {input.fa} > {output.fa} && gunzip -c {input.gtf} > {output.gtf}"

rule hisat2_index:
    input:
        f"{OUT_DIR}/00_reference/genome.fa"
    output:
        INDEX
    threads: THREADS
    conda: "env.yaml"
    shell:
        "hisat2-build -p {threads} {input} {OUT_DIR}/04_index/genome"

rule hisat2_align:
    input:
        r1 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair1.truncated.gz",
        r2 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair2.truncated.gz",
        idx = INDEX
    output:
        bam = f"{OUT_DIR}/05_alignment/{{sample}}.bam",
        bai = f"{OUT_DIR}/05_alignment/{{sample}}.bam.bai"
    threads: THREADS
    conda: "env.yaml"
    shell:
        "hisat2 -p {threads} -x {OUT_DIR}/04_index/genome -1 {input.r1} -2 {input.r2} | samtools sort -o {output.bam} && samtools index {output.bam}"
    
rule feature_counts:
    input:
        bam = expand(f"{OUT_DIR}/05_alignment/{{sample}}.bam", sample=SAMPLES),
        gtf = f"{OUT_DIR}/00_reference/annotation.gtf"
    output:
        f"{OUT_DIR}/06_counts/gene_counts.txt"
    threads: THREADS
    conda: "env.yaml"
    shell:
        "featureCounts -T {threads} -p --countReadPairs -t exon -g gene_id -a {input.gtf} -o {output} {input.bam}"