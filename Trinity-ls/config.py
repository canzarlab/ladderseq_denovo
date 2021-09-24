# Config file for denovo assembly using trinity-ls
# Authour: Shounak Chakraborty
# 23.06.2020

N_THREADS = 32


### The folder where the reads are stored.
### e.g. '/reads' (if the reads are in the folder called 'reads')
DATA_BASE_PATH = ''

# read type configurations
READS_SFX = ""
READS_TYPE = ""

# Annotation File Paths
# The path to the cdna fasta file
# e.g. 'Homo_sapiens.GRCh38.cdna.all.fa'
CDNA_FASTA_FILE = ''



# Software configurations
KALLISTO_LS_BINARY = '../ext/linuxBinaries/kallisto-ls'
SAMTOOLS_BINARY = '../ext/linuxBinaries/samtools'

# The following binaries should be in the PATH variable of your system
TRINITY_BINARY = "<PATH TO TRINITY BINARY>"



### parameters for ladder merge
BANDS = ["band1","band2","band3","band4","band5","band6","band7"]
KAL_IDX = DATA_BASE_PATH+'/kallisto_index.idx'
