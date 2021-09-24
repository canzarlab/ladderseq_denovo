from collections import defaultdict
import random

include: "config.py"
include: "Snakefile_prepareData.smk"
include: "Snakefile_runLadderMerge.smk"


rule all:
    input:
        #### denovo assembly
        expand(DATA_BASE_PATH+"/LadderSeqAssembly.fasta"),
