# SnakeMAGs
SnakeMAGs is a workflow to reconstruct prokaryotic genomes from metagenomes. The main purpose of SnakeMAGs is to process Illumina data from raw reads to metagenome-assembled genomes (MAGs).
SnakeMAGs is efficient, easy to handle and flexible to different projects. The workflow is CeCILL licensed, implemented in Snakemake (run on multiple cores) and available for Linux.
SnakeMAGs was designed to build a Termite Data Base (TDB) which contains the genomes of termite-associates micro-organisms.

![scheme of workflow](SnakeMAGs_schema.png?raw=true)

# How to use SnakeMAGs
## Install conda
The easiest way to install and run SnakeMAGs is to use [conda](https://www.anaconda.com/products/distribution). These package managers will help you to easily install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Install Snakemake
```
conda activate
conda create --prefix path/to/env/SNAKEMAKE
conda install -c bioconda snakemake --prefix path/to/env/SNAKEMAKE
conda activate path/to/env/SNAKEMAKE/
```
## Edit config file
You need to edit the config.yaml file. In particular, you need to set the correct paths and allocate the proper computational resources (threads, memory), according to your hardware. 

You also need to set the paths in all the provided .yaml files that correspond to conda environments (BEDTOOLS.yaml, COVERM.yaml, SAMTOOLS.yaml, BOWTIE2.yaml, GTDBTK.yaml, TRIMMOMATIC.yaml, BWA.yaml, IU.yaml, CHECKM.yaml, MEGAHIT.yaml, METABAT2.yaml). The user can choose these paths.

Here is an exemple of a config file:

```
#####################################################################################################
#####  _____    ___    _              _   _    ______   __    __              _______   _____   #####
##### /  ___|  |   \  | |     /\     | | / /  |  ____| |  \  /  |     /\     /  _____| /  ___|  #####
##### | (___   | |\ \ | |    /  \    | |/ /   | |____  |   \/   |    /  \    | |   __  | (___   #####
#####  \___ \  | | \ \| |   / /\ \   | |\ \   |  ____| | |\  /| |   / /\ \   | |  |_ |  \___ \  #####
#####  ____) | | |  \   |  / /__\ \  | | \ \  | |____  | | \/ | |  / /__\ \  | |____||  ____) | #####
##### |_____/  |_|   \__| /_/    \_\ |_|  \_\ |______| |_|    |_| /_/    \_\  \______/ |_____/  #####
#####                                                                                           #####
#####################################################################################################

########################
## Execution parameters
#########################

working_dir: /path/to/working/directory/                                 #The main directory for the project
raw_fastq: /usr/home/biopatic-DATA/tadrent/SnakeMAGs/raw_fastq/          #The directory that contains all the fastq files of all the samples (eg. sample1_R1.fastq & sample1_R2.fastq, sample2_R1.fastq & sample2_R2.fastq...)
suffix_1: "_1"                                                           #Main type of suffix for forword reads file (_1.fastq or _R1.fastq or _r1.fastq or _1.fq or _R1.fq or _r1.fq )
suffix_2: "_2"                                                           #Main type of suffix for reverse reads file (_2.fastq or _R2.fastq or _r2.fastq or _2.fq or _R2.fq or _r2.fq )

########################
### Conda environnemnts
##########################
IU_conda_env: "/path/to/SnakeMAGs_conda_env/IU.yaml"
TRIMMOMATIC_conda_env: "/path/to/SnakeMAGs_conda_env/TRIMMOMATIC.yaml"
BOWTIE2_conda_env: "/path/to/SnakeMAGs_conda_env/BOWTIE2.yaml"
SAMTOOLS_conda_env: "/path/to/SnakeMAGs_conda_env/SAMTOOLS.yaml"
BEDTOOLS_conda_env: "/path/to/SnakeMAGs_conda_env/BEDTOOLS.yaml"
MEGAHIT_conda_env: "/path/to/SnakeMAGs_conda_env/MEGAHIT.yaml"
BWA_conda_env: "/path/to/SnakeMAGs_conda_env/BWA.yaml"
METABAT2_conda_env: "/path/to/SnakeMAGs_conda_env/METABAT2.yaml"
CHECKM_conda_env: "/path/to/SnakeMAGs_conda_env/CHECKM.yaml"
GTDBTK_conda_env: "/path/to/SnakeMAGs_conda_env/GTDBTK.yaml"
COVERM_conda_env: "/path/to/SnakeMAGs_conda_env/COVERM.yaml"

########################
### Quality filtering
##########################
email: name.surname@your-univ.com                                        #Your e-mail addresse
threads_filter: 10                                                       #The number of threads to run this process. To be adjusted according to your hardware
ressources_filter: "mem=500, mem_mb=500000"                              #Memory according to tools need

##########################
### Adapter trimming
##########################
adapters: /path/to/working/directory/adapters.fa                         #A fasta file contanning a set of various Illumina adaptors (this file is provided and is also available on github)
trim_params: "2:40:15"                                                   #for further details, see the trimmomatic documentation
threads_trim: 10                                                         #The number of threads to run this process. To be adjusted according to your hardware
ressources_trim: "mem=500, mem_mb=500000"                                #Memory according to tools need

##########################
#### Host filtering
###########################
host_genome: "yes"                                                      #yes or no. An optional step for host-associated samples (eg. termite, human, plant...)
threads_bowtie2: 100                                                    #The number of threads to run this process. To be adjusted according to your hardware
host_genomes_directrory: /path/to/working/host_genomes/                 #the directory where the host genome is stored
host_genomes: /path/to/working/host_genomes/host_genomes.fa             #A fasta file containing the DNA sequences of the host genome(s)
threads_samtools: 100                                                   #The number of threads to run this process. To be adjusted according to your hardware
ressources_host_filtering: "mem=500, mem_mb=500000"                     #Memory according to tools need

##########################
##### Assembly
############################
threads_megahit: 150                                                   #The number of threads to run this process. To be adjusted according to your hardware
min_contig_len: 1000                                                   #Minimum length (in bp) of the assembled contigs
k_list: "21,31,41,51,61,71,81,91,99,109,119"                           #Kmer size (for further details, see the megahit documentation)
ressources_megahit: "mem=500, mem_mb=500000"                           #Memory according to tools need

##########################
###### Binning
#############################
threads_bwa: 150                                                       #The number of threads to run this process. To be adjusted according to your hardware
ressources_bwa: "mem=500, mem_mb=500000"                               #Memory according to tools need
threads_samtools: 150                                                  #The number of threads to run this process. To be adjusted according to your hardware
ressources_samtools: "mem=500, mem_mb=500000"                          #Memory according to tools need
seed: 19860615                                                         #Seed number for reproducible results
threads_metabat: 150                                                   #The number of threads to run this process. To be adjusted according to your hardware
minContig: 2500                                                        #Minimum length (in bp) of the contigs
ressources_binning: "mem=500, mem_mb=500000"                           #Memory according to tools need

##########################
####### Bins quality
##############################
threads_checkm: 150                                                    #The number of threads to run this process. To be adjusted according to your hardware
ressources_checkm: "mem=500, mem_mb=500000"                            #Memory according to tools need

##########################
######## Classification
###############################
threads_gtdb: 16                                                       #The number of threads to run this process. To be adjusted according to your hardware
ressources_gtdb: "mem=500, mem_mb=500000"                              #Memory according to tools need

##########################
####### Abundances
##############################
threads_coverM: 16                                                     #The number of threads to run this process. To be adjusted according to your hardware
ressources_coverM: 500                                                 #Memory according to tools need
```
