import os
configfile: "config.yaml"

OUTPUT_DIR = config["OUTPUT_DIR"]
ANALYSIS_DIR = config["ANALYSIS_DIR"]
PHENOTYPES = config["phenotype"]
GENE_BED_DIR = config["GENE_BED_DIR"]
GENE_BEDS = config["gene_beds"]
SCRIPT = config["script"]
SAMPLE_TO_FETCH_READS_FROM = os.path.basename(config["SAMPLE_TO_FETCH_READS_FROM"])

settings = f"mac{config['k_min_allowed_appearance_count']}p{int(config['k_min_percent_appearance_in_each_strand_form']*100)}"

# extensions
ext_kmers_to_use = ["", ".no_pass_kmers", ".shareness", ".stats.both", ".stats.only_canonical", ".stats.only_non_canonical"]
ext_kmers_table = [".table", ".names"]
ext_bwa_ref_index = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]


rule all:
    input:
        expand(f"{ANALYSIS_DIR}/{{phenotype}}/reads/reads.intersect.{{gene_beds}}", phenotype=PHENOTYPES, gene_beds=GENE_BEDS),
        expand(f"{ANALYSIS_DIR}/{{phenotype}}/reads/reads.closest.{{gene_beds}}", phenotype=PHENOTYPES, gene_beds=GENE_BEDS),
        expand(f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.intersect.{{gene_beds}}", phenotype=PHENOTYPES, gene_beds=GENE_BEDS),
        expand(f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.closest.{{gene_beds}}", phenotype=PHENOTYPES, gene_beds=GENE_BEDS),



# prepare phenotype value file and kmers_list_paths.txt for downstream kmerGWAS
rule prepare_required_files:
    output:
        pheno_file = f"{OUTPUT_DIR}/{{phenotype}}/{{phenotype}}_phenotype_values.tsv",
        kmers_list_paths_file = f"{OUTPUT_DIR}/{{phenotype}}/kmers_list_paths.txt"
    params:
        outdir = f"{OUTPUT_DIR}/{{phenotype}}",
        script = SCRIPT,
        phenotype = "{phenotype}",
        extra = "-rf"
    conda:
        "envs/kmerGWAS.yaml"
    shell:
        "python {params.script} PrepFiles --xlsx {config[phenotype_book]} --phenotype {params.phenotype} --sample_csv {config[sample_csv]} --strand_dir {config[STRAND_DIR]} -o {params.outdir} {params.extra}"

# >>> kmerGWAS <<<

rule list_kmers_found_in_multiple_samples:
    input:
        f"{OUTPUT_DIR}/{{phenotype}}/kmers_list_paths.txt"
    params:
        k = config["k_length"],
        mac = config["k_min_allowed_appearance_count"],
        p = config["k_min_percent_appearance_in_each_strand_form"],
        kmers_to_use_prefix = f"{OUTPUT_DIR}/{{phenotype}}/{settings}.kmers_to_use"
    output:
        expand(f"{OUTPUT_DIR}/"+"{{phenotype}}"+f"/{settings}.kmers_to_use"+"{ext}", ext=ext_kmers_to_use)
    conda:
        "envs/kmerGWAS.yaml"
    shell:
        "{config[LIB_DIR]}/bin/list_kmers_found_in_multiple_samples \
        -l {input} \
        -k {params.k} \
        --mac {params.mac} \
        -p {params.p} \
        -o {params.kmers_to_use_prefix}"


rule create_kmers_table:
    input:
        kmers_list_paths_file = f"{OUTPUT_DIR}/{{phenotype}}/kmers_list_paths.txt",
        kmers_to_use = f"{OUTPUT_DIR}/{{phenotype}}/{settings}.kmers_to_use"
    params:
        k = config["k_length"],
        kmers_table_prefix = f"{OUTPUT_DIR}/{{phenotype}}/{settings}.kmers_table"
    output:
        expand(f"{OUTPUT_DIR}/"+"{{phenotype}}"+f"/{settings}.kmers_table"+"{ext}", ext=ext_kmers_table)
    conda:
        "envs/kmerGWAS.yaml"
    shell:
        "{config[LIB_DIR]}/bin/build_kmers_table \
        -l {input.kmers_list_paths_file} \
        -k {params.k} \
        -a {input.kmers_to_use} \
        -o {params.kmers_table_prefix}"


rule calculate_kinship_matrix:
    input:
        expand(f"{OUTPUT_DIR}/"+"{{phenotype}}"+f"/{settings}.kmers_table"+"{ext}", ext=ext_kmers_table)
    params:
        k = config["k_length"],
        maf = config["k_min_allowed_appearance_frequency"],
        kmers_table_prefix = f"{OUTPUT_DIR}/{{phenotype}}/{settings}.kmers_table"
    output:
        f"{OUTPUT_DIR}/{{phenotype}}/{settings}.kmers_table.kinship"
    conda:
        "envs/kmerGWAS.yaml"
    shell:
        "{config[LIB_DIR]}/bin/emma_kinship_kmers \
        -t {params.kmers_table_prefix} \
        -k {params.k} \
        --maf {params.maf} > {output}"
    
    

rule kmerGWAS:
    input:
        pheno_file = f"{OUTPUT_DIR}/{{phenotype}}/{{phenotype}}_phenotype_values.tsv",
        kmers_table = expand(f"{OUTPUT_DIR}/"+"{{phenotype}}"+f"/{settings}.kmers_table"+"{ext}", ext=ext_kmers_table),
        kinship = f"{OUTPUT_DIR}/{{phenotype}}/{settings}.kmers_table.kinship"
    params:
        k = config["k_length"],
        kmers_table_prefix = f"{OUTPUT_DIR}/{{phenotype}}/{settings}.kmers_table",
        results_folder = f"{OUTPUT_DIR}/{{phenotype}}/results"
    output:
        f"{OUTPUT_DIR}/{{phenotype}}/results",
        f"{OUTPUT_DIR}/{{phenotype}}/results/kmers/pass_threshold_5per"
    threads:
        config["THREADS"]
    conda:
        "envs/kmerGWAS.yaml"
    shell:
        "rm -r {params.results_folder} && "
        "export LD_LIBRARY_PATH=/g/data/xf3/miniconda/envs/kmerGWAS/lib && "
        "python2.7 {config[LIB_DIR]}/kmers_gwas.py \
        --pheno {input.pheno_file} \
        --kmers_table {params.kmers_table_prefix} \
        -l {params.k} \
        -p {threads} \
        --outdir {output[0]}"


rule kmers2fasta:
    input:
        f"{OUTPUT_DIR}/{{phenotype}}/results/kmers/pass_threshold_5per"
    params:
        script = SCRIPT,
        phenotype = "{phenotype}",
        outdir = f"{ANALYSIS_DIR}/{{phenotype}}"
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/{{phenotype}}_kmers.fasta"
    conda:
        "envs/kmerGWAS.yaml"
    shell:
        "python {params.script} kmers2fasta --kmers_in {input} --phenotype {params.phenotype} -o {params.outdir}"


rule fetch_reads_with_kmers:
    input:
        kmers_fasta = f"{ANALYSIS_DIR}/{{phenotype}}/{{phenotype}}_kmers.fasta",
        r1 = f"{config['SAMPLE_TO_FETCH_READS_FROM']}_1.fastq.gz",
        r2 = f"{config['SAMPLE_TO_FETCH_READS_FROM']}_2.fastq.gz"
    params:
        k = config["k_length"],
        outbase = f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}"
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}_R1.fastq",
        f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}_R2.fastq"
    shell:
        "{config[FETCH_READS_WITH_KMERS_DIR]}/fetch_reads {input.r1} {input.r2} {input.kmers_fasta} {params.k} {params.outbase}"



rule bwamem2_index_genome:
    input:
        config["ref_genome"]
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "bwa-mem2 index -p {input}"


rule bedtools_sort_gene_beds:
    input:
        f"{GENE_BED_DIR}/{{gene_beds}}"
    output:
        f"{GENE_BED_DIR}/sorted.{{gene_beds}}"
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "bedtools sort -i {input} > {output}"



### READS MAPPING ###


rule READS_bwamem2_map:
    input:
        r1 = f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}_R1.fastq",
        r2 = f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}_R2.fastq"
    params:
        config["ref_genome"]
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}.{config['ref_genome_name']}.bam"
    conda:
        "envs/kmerAnalysis.yaml"
    threads:
        config["THREADS"]
    shell:
        "bwa-mem2 mem -t {threads} {params} {input.r1} {input.r2} | samtools sort -O BAM -@ {threads} -o {output} - "


rule READS_keep_only_mapped_reads:
    input:
        f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}.{config['ref_genome_name']}.bam"
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}.{config['ref_genome_name']}.mapped_only.bam"
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "samtools view -b -F 4 {input} > {output}"


rule READS_bedtools_bamtobed:
    input:
        f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}.{config['ref_genome_name']}.mapped_only.bam"
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}.{config['ref_genome_name']}.mapped_only.bam.bed"
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "bedtools bamtobed -i {input} | bedtools sort > {output}"

        

rule READS_bedtools_intersect:
    input:
        readbed = f"{ANALYSIS_DIR}/{{phenotype}}/reads/fetched.{SAMPLE_TO_FETCH_READS_FROM}.{config['ref_genome_name']}.mapped_only.bam.bed",
        genebed = f"{GENE_BED_DIR}/sorted.{{gene_beds}}"
    output:
        intersect = f"{ANALYSIS_DIR}/{{phenotype}}/reads/reads.intersect.{{gene_beds}}",
        closest = f"{ANALYSIS_DIR}/{{phenotype}}/reads/reads.closest.{{gene_beds}}"
    params:
        tmp = f"{ANALYSIS_DIR}/{{phenotype}}/reads/tmp"
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "bedtools intersect -a {input.readbed} -b {input.genebed} | sort > {output.intersect} && "
        "bedtools closest -a {input.readbed} -b {input.genebed} -io | sort > {params.tmp} && "
        "cut {params.tmp} -f7,8,9,10,11,12 | sort | uniq | tail -n +2 > {output.closest} && "
        "rm {params.tmp}"


### KMERS MAPPING ###


rule KMERS_bwamem2_map:
    input:
        f"{ANALYSIS_DIR}/{{phenotype}}/{{phenotype}}_kmers.fasta"
    params:
        config["ref_genome"]
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.{config['ref_genome_name']}.bam"
    conda:
        "envs/kmerAnalysis.yaml"
    threads:
        config["THREADS"]
    shell:
        "bwa-mem2 mem -t {threads} {params} {input} | samtools sort -O BAM -@ {threads} -o {output} - "


rule KMERS_keep_only_mapped_kmers:
    input:
        f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.{config['ref_genome_name']}.bam"
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.{config['ref_genome_name']}.mapped_only.bam"
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "samtools view -b -F 4 {input} > {output}"


rule KMERS_bedtools_bamtobed:
    input:
        f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.{config['ref_genome_name']}.mapped_only.bam"
    output:
        f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.{config['ref_genome_name']}.mapped_only.bam.bed"
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "bedtools bamtobed -i {input} | bedtools sort > {output}"


rule KMERS_bedtools_intersect:
    input:
        kmerbed = f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.{config['ref_genome_name']}.mapped_only.bam.bed",
        genebed = f"{GENE_BED_DIR}/sorted.{{gene_beds}}"
    output:
        intersect = f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.intersect.{{gene_beds}}",
        closest = f"{ANALYSIS_DIR}/{{phenotype}}/kmers/kmers.closest.{{gene_beds}}"
    params:
        tmp = f"{ANALYSIS_DIR}/{{phenotype}}/reads/tmp"
    conda:
        "envs/kmerAnalysis.yaml"
    shell:
        "bedtools intersect -a {input.kmerbed} -b {input.genebed} | sort > {output.intersect} && "
        "bedtools closest -a {input.kmerbed} -b {input.genebed} -io | sort > {params.tmp} && "
        "cut {params.tmp} -f7,8,9,10,11,12 | sort | uniq | tail -n +2 > {output.closest} && "
        "rm {params.tmp}"