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
        f"{OUT_DIR}/06_counts/gene_counts.txt",
        expand(f"{OUT_DIR}/07_kallisto/{{sample}}/abundance.tsv", sample=SAMPLES),
        expand(f"{OUT_DIR}/08_rseqc/{{sample}}.bam_stat.txt", sample=SAMPLES),
        expand(f"{OUT_DIR}/08_rseqc/{{sample}}.infer_experiment.txt", sample=SAMPLES),
        expand(f"{OUT_DIR}/08_rseqc/{{sample}}.read_distribution.txt", sample=SAMPLES),
        expand(f"{OUT_DIR}/08_rseqc/{{sample}}.inner_distance_freq.txt", sample=SAMPLES),
        f"{OUT_DIR}/08_rseqc/all_samples.geneBodyCoverage.curves.pdf",
        f"{OUT_DIR}/09_deeptools/correlation_heatmap.pdf",
        f"{OUT_DIR}/09_deeptools/pca_plot.pdf",


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
        gtf = f"{OUT_DIR}/00_reference/annotation.gtf",
        bed = f"{OUT_DIR}/00_reference/annotation.bed"
    conda: "env.yaml"
    shell:
        "gunzip -c {input.fa} > {output.fa} && gunzip -c {input.gtf} > {output.gtf} &&"
        "gtfToGenePred {output.gtf} {OUT_DIR}/00_reference/temp.gp && genePredToBed {OUT_DIR}/00_reference/temp.gp {output.bed} && rm {OUT_DIR}/00_reference/temp.gp"

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

#--Alignment-free quantification-----------------------------------------

rule kallisto_index:
    input:
        config["paths"]["transcriptome_gz"]
    output:
        f"{OUT_DIR}/04_index/transcriptome.idx"
    conda: "env.yaml"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto_quant:
    input:
        r1 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair1.truncated.gz",
        r2 = f"{OUT_DIR}/02_trimmed/{{sample}}.pair2.truncated.gz",
        idx = f"{OUT_DIR}/04_index/transcriptome.idx"
    output:
        f"{OUT_DIR}/07_kallisto/{{sample}}/abundance.tsv"
    threads: THREADS
    conda: "env.yaml"
    shell:
        "kallisto quant -i {input.idx} -o {OUT_DIR}/07_kallisto/{wildcards.sample} -t {threads} {input.r1} {input.r2}"

#--ALIGNMENT QUALITY ASSESSMENT-------------------------------------

rule rseqc_sample:
    input:
        bam = f"{OUT_DIR}/05_alignment/{{sample}}.bam",
        bai = f"{OUT_DIR}/05_alignment/{{sample}}.bam.bai",
        bed = f"{OUT_DIR}/00_reference/annotation.bed"
    output:
        stat = f"{OUT_DIR}/08_rseqc/{{sample}}.bam_stat.txt",
        strand = f"{OUT_DIR}/08_rseqc/{{sample}}.infer_experiment.txt",
        dist = f"{OUT_DIR}/08_rseqc/{{sample}}.read_distribution.txt",
        inner = f"{OUT_DIR}/08_rseqc/{{sample}}.inner_distance_freq.txt"
    conda: "env.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output.stat} && "
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output.strand} && read_distribution.py -r {input.bed} -i {input.bam} > {output.dist} && "
        "inner_distance.py -r {input.bed} -i {input.bam} -o {OUT_DIR}/08_rseqc/{wildcards.sample}"

rule rseqc_genebody:
    input:
        bams = expand(f"{OUT_DIR}/05_alignment/{{sample}}.bam", sample=SAMPLES),
        bais = expand(f"{OUT_DIR}/05_alignment/{{sample}}.bam.bai", sample=SAMPLES),
        bed = f"{OUT_DIR}/00_reference/annotation.bed"
    output:
        f"{OUT_DIR}/08_rseqc/all_samples.geneBodyCoverage.curves.pdf"
    params:
        bam_csv = ",".join(expand(f"{OUT_DIR}/05_alignment/{{sample}}.bam", sample=SAMPLES))
    conda: "env.yaml"
    shell:
        "geneBody_coverage.py -r {input.bed} -i {params.bam_csv} -o {OUT_DIR}/08_rseqc/all_samples"

rule deeptools_qc:
    input:
        bams = expand(f"{OUT_DIR}/05_alignment/{{sample}}.bam", sample=SAMPLES),
        bais = expand(f"{OUT_DIR}/05_alignment/{{sample}}.bam.bai", sample=SAMPLES)
    output:
        npz = f"{OUT_DIR}/09_deeptools/results.npz",
        corr = f"{OUT_DIR}/09_deeptools/correlation_heatmap.pdf",
        pca = f"{OUT_DIR}/09_deeptools/pca_plot.pdf"
    threads: THREADS
    conda: "env.yaml"
    shell:
        "multiBamSummary bins --bamfiles {input.bams} -p {threads} -o {output.npz} && "
        "plotCorrelation -in {output.npz} --corMethod pearson --whatToPlot heatmap --plotNumbers -o {output.corr} && "
        "plotPCA -in {output.npz} -o {output.pca}"