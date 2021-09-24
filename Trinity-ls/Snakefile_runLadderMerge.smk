### Snakemake file for ###
### 1. Denovo assembly using trinity ###

### Author: Shounak ###
### 23.06.2020 ###

# FILTER_EXTENT = "CORRECT_PLUS_ONE"
# FILTER_RES = "HIGH"
# MOD_KALLISTO_FILTER = 0.1


rule ladder_seq_trinity_individual:
    input:
        DATA_BASE_PATH+"/{band}/R1" + READS_SFX + "." + READS_TYPE,
        DATA_BASE_PATH+"/{band}/R2" + READS_SFX + "." + READS_TYPE
    threads:
    	N_THREADS
    output:
        "{DATA_BASE_PATH}/{band}_trinity/Trinity.fasta"
    message: "Running Trinity on {wildcards.band} ..."
    shell:
        """
        mkdir -p {DATA_BASE_PATH}/{wildcards.band}_trinity/
        {TRINITY_BINARY} --seqType {READS_TYPE} --left {input[0]} --right {input[1]} --CPU 16 --max_memory 20G --output {DATA_BASE_PATH}/{wildcards.band}_trinity
        """

rule filter_trinity_individual:
    input:
        DATA_BASE_PATH+"/{band}_trinity/Trinity.fasta",
        DATA_BASE_PATH+"/migProb.tsv"
    threads:
    	N_THREADS
    output:
        "{DATA_BASE_PATH}/{band}/Trinity_filtered.fasta"
    message: "Filtering Trinity output on {wildcards.band}..."
    shell:
        """
        python FilterLadderFasta.py {input[0]} {output[0]} length {wildcards.band} CORRECT_PLUS_ONE HIGH {input[1]}
        """


rule concat_contigs:
    input:
        expand(DATA_BASE_PATH+"/{band}/"+"Trinity_filtered.fasta", band=BANDS)
    threads:
        N_THREADS
    output:
        "{DATA_BASE_PATH}/concat/concatenated.fasta"
    shell:
        """
        mkdir -p {DATA_BASE_PATH}/concat
        cat {input} | grep -v '^#' > {DATA_BASE_PATH}/concat/concatenated.fasta
        """


rule run_kallistoMod_index:
	input:
		DATA_BASE_PATH+"/concat/concatenated.fasta"
	threads:
	        N_THREADS
	output:
		"{DATA_BASE_PATH}/ladderMerge/modKallistoIndex.idx"
	shell:
		"""
		mkdir -p {DATA_BASE_PATH}/ladderMerge
		{KALLISTO_LS_BINARY} index --make-unique -i {output[0]} {input[0]}
		"""


rule run_kallisto_modified_denovo:
	input:
		kalIndex = DATA_BASE_PATH+"/ladderMerge/modKallistoIndex.idx",
		ladderReads = expand(DATA_BASE_PATH+"/{band}/R{pair}_sim.fq",band=BANDS,pair=[1,2])
	threads:
	        N_THREADS
	output:
		"{DATA_BASE_PATH}/ladderMerge/abundance.tsv"
	shell:
		"""
		{KALLISTO_LS_BINARY} quant -t {threads} -o {DATA_BASE_PATH}/ladderMerge/ -i {input.kalIndex} {input.ladderReads}
		"""


rule run_filterModKallistoGff:
    input:
        DATA_BASE_PATH+"/concat/concatenated.fasta",
        DATA_BASE_PATH+"/ladderMerge/abundance.tsv"
    threads:
        N_THREADS
    output:
        "{DATA_BASE_PATH}/ladderMerge/transcripts_tobeKept.txt",
        "{DATA_BASE_PATH}/ladderMerge/LadderSeqAssembly.fasta"
    shell:
    ## the fourth and 5th argument is filterByTpm and filterByCount
        """
        awk '$5 >= 0.1 {{print $1}}' {input[1]} > {output[0]}
        python FilterLadderFasta.py {input[0]} {output[1]} id {output[0]}
        """


rule cleanup:
	input:
		DATA_BASE_PATH+"/ladderMerge/LadderSeqAssembly.fasta"
	threads:
	        N_THREADS
	output:
		DATA_BASE_PATH+"/LadderSeqAssembly.fasta"
	shell:
		"""
        mv {input} {DATA_BASE_PATH}/
        rm -r {DATA_BASE_PATH}/ladderMerge
        rm -r {DATA_BASE_PATH}/band*
        rm -r {DATA_BASE_PATH}/original
        rm -r {DATA_BASE_PATH}/merged
		"""
