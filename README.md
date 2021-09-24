fkallisto# *De novo* assembly from Ladder-seq

![Ladder-Seq Assembly pipeline](/AssemblyImage400.png)


**Ladder-seq**  is the concerted advancement of the RNA-seq protocol and its computational methods. It experimentally separates transcripts according to their length prior to sequencing to achieve a "coloring" of reads that connects them along transcript isoforms.

In the ***de novo* transcript assembly** from Ladder-seq reads we use [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) to build a separate de Bruijn graph in each band and apply transcript length constraints derived from read colorings to break too long erroneous fusions and to eliminate too short transcript fragments. Finally, reads are assigned to assembled transcripts guided by their coloring using [kallisto-ls](https://github.com/canzarlab/ladderseq_quant), our extension of [kallisto](https://pachterlab.github.io/kallisto/about) to Ladder-seq.

 <br />

## Installation

This repository includes precompiled [samtools 1.10](https://github.com/samtools/samtools/releases) and kallisto-ls binaries for Linux x86_64.

The Trinity-ls pipeline uses the de novo transcriptome assembler Trinity which can be built from source as described [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Installing-Trinity). The path to the Trinity binary needs to be adjusted in the `config.py` file (see below).

Similarly, kallisto-ls can be built from source as described [here](https://github.com/canzarlab/ladderseq_quant) and the path to the kallisto-ls binary can be adjusted in the `config.py` file (see below).

All required steps, including the estimation of migration patterns, length-contrained transcript assembly, integration, and quantification are run by the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system, which needs to be [installed](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on your system.

## *De novo* transcript assembly using Trinity

We refer to the Trinity based assembly workflow for Ladder-seq as __Trinity-ls__. It assumes that Ladder-seq input reads for each band are stored in a common directory with the following file name convention: ```R1_band<i>.fq``` and ```R2_band<i>.fq``` for paired-end reads in band ```i```. __Trinity-ls__ currently assumes a separation of reads into 7 bands.

In file `config.py` in the Trinity-ls subdirectory, the user needs to set variables for the directory path to the input read files and to cDNA sequences of spike-in RNA transcripts or supplemental cell line transcripts in FASTA format. More detailed instructions are provided in file `config.py`. To run __Trinity-ls__, simply execute the Snakemake workflow:

```shell
cd Trinity-ls
snakemake -s Snakemake_runPipeline.smk
```
The number of cores to be used by Snakemake can be specified using `-j <number of cores>`, which is described in more detailed in the Snakemake [documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options).

Assembled transcripts are reported in file `LadderSeqAssembly.gtf` in the directory specified by variable `DATA_BASE_PATH` in `config.py`.

## Transcript assembly using an alternative *de novo* assembler

In directory `GenericAssembler` we provide a more generic workflow in which Trinity can easily be replaced by a *de novo* RNA-seq assembly method of the user's choice.
In file `Snakefile_runLadderMerge.smk` the user needs to replace a single occurrences of
```shell
<YOUR CODE GOES HERE>
```
with the command used to call the preferred assemby method, using `{input[0]}` and `{input[1]}` to specify its input paired-end read files in `.fastq` or `.fasta` format, and `{output[0]}` to specify the output file containing transcripts in `.gtf` format:

As described above, set variables in `config.py` in directory `GenericAssembler` and run the Snakemake workflow using

```shell
cd GenericAssembler
snakemake -s Snakemake_runPipeline.smk
```
