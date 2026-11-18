configfile: "config.yaml"

RAW_DIR = config["paths"]["raw_dir"]
OUT_DIR = config["paths"]["out_dir"]
SAMPLES = glob_wildcards(f"{RAW_DIR}/{{sample}}_1.fastq").sample

rule all:
    input:
        expand(f"{OUT_DIR}/01_qc_raw/{{sample}}_{{read}}_fastqc.html",sample=SAMPLES, read=[1, 2]),
        expand(f"{OUT_DIR}/03_qc_trimmed/{{sample}}.pair{{read}}.truncated_fastqc.html", sample=SAMPLES, read=[1, 2])

rule fastqc_raw:
    input:
        r1 = f"{RAW_DIR}/{{sample}}_1.fastq",
        r2 = f"{RAW_DIR}/{{sample}}_2.fastq"
    output:
        f"{OUT_DIR}/01_qc_raw/{{sample}}_1_fastqc.html",
        f"{OUT_DIR}/01_qc_raw/{{sample}}_2_fastqc.html"
    threads: config["params"]["threads"]
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
    threads: config["params"]["threads"]
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
    threads: config["params"]["threads"]
    conda: "env.yaml"
    shell:
        "fastqc -t {threads} -o {OUT_DIR}/03_qc_trimmed {input.r1} {input.r2}"
