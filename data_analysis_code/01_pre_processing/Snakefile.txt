configfile: "config.yaml"

rule all:
    input:
        expand(expand("scratch/combined/{plate}_{primer}.filt.fwd.fq.gz", plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"]),
        expand(expand("scratch/combined/{plate}_{primer}.filt.rev.fq.gz", plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"]),
        "scratch/combined/ASVs.msa.treefile",
        "scratch/combined/phyloseq_taxa.rds"

rule demultiplex:
    conda:
        "cutadapt"
    input:
        golay_barcodes="data/golay_linker_barcodes_combined.fasta",
        golay_barcodes_reverse_order="data/golay_linker_barcodes_reverse_combined.fasta",
        fwd="data/{plate}_R1_001.fastq.gz",
        rev="data/{plate}_R2_001.fastq.gz"
    output:
        expand("scratch/{plate}/{primer}.fwd.fq.gz", zip, primer=config["PRIMERS"], allow_missing=True),
        expand("scratch/{plate}/{primer}.rev.fq.gz", zip, primer=config["PRIMERS"], allow_missing=True)
    params:
        primer_combos=lambda wildcards:'{name}',
        folder=directory("scratch/{plate}")
    threads: 1
    log:
        "logs/{plate}_demultiplex_reads.log"
    benchmark:
        "benchmarks/{plate}_demultiplex_reads.tsv"
    shell:
        """
        mkdir -p {params.folder}
        cutadapt -j {threads} -e 2 -m 1 --pair-adapters -g ^file:{input.golay_barcodes} -G ^file:{input.golay_barcodes_reverse_order} -o {params.folder}/{params.primer_combos}.fwd.fq.gz -p {params.folder}/{params.primer_combos}.rev.fq.gz {input.fwd} {input.rev} 2>&1 | tee {log}
        """
rule FilterAndTrim:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        fwd="scratch/{plate}/{primer}.fwd.fq.gz",
        rev="scratch/{plate}/{primer}.rev.fq.gz"
    output:
        fwd="scratch/combined/{plate}_{primer}.filt.fwd.fq.gz",
        rev="scratch/combined/{plate}_{primer}.filt.rev.fq.gz"
    threads: 1
    log:
        "logs/{plate}/{primer}_demultiplex_reads.log"
    benchmark:
        "benchmarks/{plate}/{primer}_demultiplex_reads.tsv"
    shell:
        """
        Rscript scripts/FilterAndTrim_func.R {input.fwd} {input.rev} {output.fwd} {output.rev} 2>&1 | tee {log}
        touch {output.fwd} {output.rev} # These files are not made if input is poor, but snakemake will quit if it doesn't get an output
        """
rule EstimateErrorsFwd1:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errF_1.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsFwd1.log"
    benchmark:
        "benchmarks/EstimateErrorsFwd1.tsv"
    shell:
        """
        Rscript scripts/EC1_fwd.R {params} 2>&1 | tee {log}
        """

rule EstimateErrorsRev1:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errR_1.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsRev1.log"
    benchmark:
        "benchmarks/EstimateErrorsRev1.tsv"
    shell:
        """
        Rscript scripts/EC1_rev.R {params} 2>&1 | tee {log}
        """
rule EstimateErrorsFwd2:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errF_2.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsFwd2.log"
    benchmark:
        "benchmarks/EstimateErrorsFwd2.tsv"
    shell:
        """
        Rscript scripts/EC2_fwd.R {params} 2>&1 | tee {log}
        """

rule EstimateErrorsRev2:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errR_2.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsRev2.log"
    benchmark:
        "benchmarks/EstimateErrorsRev2.tsv"
    shell:
        """
        Rscript scripts/EC2_rev.R {params} 2>&1 | tee {log}
        """

rule EstimateErrorsFwd3:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errF_3.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsFwd3.log"
    benchmark:
        "benchmarks/EstimateErrorsFwd3.tsv"
    shell:
        """
        Rscript scripts/EC3_fwd.R {params} 2>&1 | tee {log}
        """

rule EstimateErrorsRev3:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errR_3.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsRev3.log"
    benchmark:
        "benchmarks/EstimateErrorsRev3.tsv"
    shell:
        """
        Rscript scripts/EC3_rev.R {params} 2>&1 | tee {log}
        """

rule EstimateErrorsFwd4:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errF_4.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsFwd4.log"
    benchmark:
        "benchmarks/EstimateErrorsFwd4.tsv"
    shell:
        """
        Rscript scripts/EC4_fwd.R {params} 2>&1 | tee {log}
        """

rule EstimateErrorsRev4:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"])
    output:
        "scratch/combined/errR_4.rds"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/EstimateErrorsRev4.log"
    benchmark:
        "benchmarks/EstimateErrorsRev4.tsv"
    shell:
        """
        Rscript scripts/EC4_rev.R {params} 2>&1 | tee {log}
        """
rule dada2:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        expand(expand(rules.FilterAndTrim.output, plate=config["PLATES"], allow_missing=True), zip, primer=config["PRIMERS"]),
        rules.EstimateErrorsFwd1.output,
        rules.EstimateErrorsRev1.output,
        rules.EstimateErrorsFwd2.output,
        rules.EstimateErrorsRev2.output,
        rules.EstimateErrorsFwd3.output,
        rules.EstimateErrorsRev3.output,
        fwd4=rules.EstimateErrorsFwd4.output,
        rev4=rules.EstimateErrorsRev4.output,
        primer_combos="data/primer_combos.csv"
    output:
        "scratch/combined/phyloseq.rds",
        "scratch/combined/ASVs.fna"
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/dada2.log"
    benchmark:
        "benchmarks/dada2.tsv"
    shell:
        """
        Rscript scripts/dada2.R {params} {input.fwd4} {input.rev4} {input.primer_combos} 2>&1 | tee {log}
        """

rule BuildTree:
    conda:
        "iqtree"
    input:
        "scratch/combined/ASVs.fna"
    output:
        msa="scratch/combined/ASVs.msa",
        tree="scratch/combined/ASVs.msa.treefile"
    params:
        folder=directory("scratch/combined"),
        msa="ASVs.msa"
    threads: 16
    resources:
        mem_gb=120
    log:
        "scratch/combined/ASVs.msa.log"
    benchmark:
        "benchmarks/build_tree.tsv"
    shell:
        """
        mafft --auto --thread {threads} {input} > {output.msa}
        cd {params.folder}
        iqtree -s {params.msa} --seed 1234 -T AUTO --mem {resources.mem_gb}G -m MFP
        """

rule AssignTaxonomy:
    conda:
        "R"
        #"envs/dada2.yaml"
    input:
        ps="scratch/combined/phyloseq.rds",
        silva="data/silva_nr99_v138.1_train_set.fa.gz",
        silva_species="data/silva_species_assignment_v138.1.fa.gz"
    output:
        "scratch/combined/phyloseq_taxa.rds",
    params:
        directory("scratch/combined")
    threads: 16
    log:
        "logs/AssignTaxonomy.log"
    benchmark:
        "benchmarks/AssignTaxonomy.tsv"
    shell:
        """
        Rscript scripts/assignTaxonomy.R {params} {input.ps} {input.silva} {input.silva_species} 2>&1 | tee {log}
        """
