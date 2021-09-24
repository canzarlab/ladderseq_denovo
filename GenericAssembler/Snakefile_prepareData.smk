
rule kallisto_index:
    input:
        CDNA_FASTA_FILE
    output:
        KAL_IDX
    threads: 1
    shell:
        KALLISTO_LS_BINARY + ' index '
        '-i ' + KAL_IDX + ' ' +
        CDNA_FASTA_FILE

rule organizeReads:
    input:
        DATA_BASE_PATH+"/R1_{individualBand}.fq",
        DATA_BASE_PATH+"/R2_{individualBand}.fq"
    threads:
        N_THREADS
    output:
        DATA_BASE_PATH+"/{individualBand}/R1_sim.fq",
        DATA_BASE_PATH+"/{individualBand}/R2_sim.fq"
    run:
        shell("mkdir -p {DATA_BASE_PATH}/{wildcards.individualBand}/")
        shell("cp {input[0]} {DATA_BASE_PATH}/{wildcards.individualBand}/R1_sim.fq")
        shell("cp {input[1]} {DATA_BASE_PATH}/{wildcards.individualBand}/R2_sim.fq")


rule estimateMigProbs:
    input:
        KAL_IDX,
        ladderReads = expand(DATA_BASE_PATH+"/{band}/R{pair}_sim.fq",band=BANDS,pair=[1,2])
    threads:
        N_THREADS
    output:
        DATA_BASE_PATH+"/migProb.tsv"
    run:
        shell("{KALLISTO_LS_BINARY} migrate -o {DATA_BASE_PATH}/ -i {input}")
