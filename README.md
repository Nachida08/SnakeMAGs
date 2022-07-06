# About SnakeMAGs
SnakeMAGs is a workflow to reconstruct prokaryotic genomes from metagenomes. The main purpose of SnakeMAGs is to process Illumina data from raw reads to metagenome-assembled genomes (MAGs).
SnakeMAGs is efficient, easy to handle and flexible to different projects. The workflow is CeCILL licensed, implemented in Snakemake (run on multiple cores) and available for Linux.

![scheme of workflow](SnakeMAGs_schema.jpg?raw=true)

# How to use SnakeMAGs
## Install conda
The easiest way to install and run SnakeMAGs is to use [conda](https://www.anaconda.com/products/distribution). These package managers will help you to easily install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Install and activate Snakemake environement
```
conda activate
conda create --prefix path/to/env/SNAKEMAKE
conda install -c bioconda snakemake --prefix path/to/env/SNAKEMAKE
conda activate path/to/env/SNAKEMAKE/
```
## SnakeMAGs input files
- Illumina paired-end reads in FASTQ.
- Adapter sequence file ([adapter.fa](https://github.com/Nachida08/SnakeMAGs/blob/main/adapters.fa)).
- Host genome sequences in FASTA (if host_genome: "yes")

## Download GENOME TAXONOMY DATABASE (GTDB)
GTDB-Tk requires ~66G+ of external data (GTDB) that need to be downloaded and unarchived.
```
#Download the latest release (tested with release207)
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
#Decompress
tar -xzvf *tar.gz
```
All you have to do now is to indicate the path to the data base folder in the config file, Classification section.

## Edit config file
You need to edit the config.yaml file. In particular, you need to set the correct paths and allocate the proper computational resources (threads, memory), according to your hardware. 

You also need to set the paths in all the provided .yaml files (available in [SnakeMAGs_conda_env directory](https://github.com/Nachida08/SnakeMAGs/tree/main/SnakeMAGs_conda_env)) that correspond to conda environments (BEDTOOLS.yaml, COVERM.yaml, SAMTOOLS.yaml, BOWTIE2.yaml, GTDBTK.yaml, TRIMMOMATIC.yaml, BWA.yaml, IU.yaml, CHECKM.yaml, MEGAHIT.yaml, METABAT2.yaml). The user can choose these paths.

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

############################
### Execution parameters ###
############################

working_dir: /path/to/working/directory/                                 #The main directory for the project
raw_fastq: /path/to/raw_fastq/                                           #The directory that contains all the fastq files of all the samples (eg. sample1_R1.fastq & sample1_R2.fastq, sample2_R1.fastq & sample2_R2.fastq...)
suffix_1: "_1"                                                           #Main type of suffix for forword reads file (_1.fastq or _R1.fastq or _r1.fastq or _1.fq or _R1.fq or _r1.fq )
suffix_2: "_2"                                                           #Main type of suffix for reverse reads file (_2.fastq or _R2.fastq or _r2.fastq or _2.fq or _R2.fq or _r2.fq )

###########################
### Conda environnemnts ###
###########################
 
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

#########################
### Quality filtering ###
#########################
email: name.surname@your-univ.com                                        #Your e-mail address
threads_filter: 10                                                       #The number of threads to run this process. To be adjusted according to your hardware
ressources_filter: 500                                                   #Memory according to tools need

########################
### Adapter trimming ###
########################
adapters: /path/to/working/directory/adapters.fa                         #A fasta file contanning a set of various Illumina adaptors (this file is provided and is also available on github)
trim_params: "2:40:15"                                                   #For further details, see the trimmomatic documentation
threads_trim: 10                                                         #The number of threads to run this process. To be adjusted according to your hardware
ressources_trim: 500                                #Memory according to tools need

######################
### Host filtering ###
######################
host_genome: "yes"                                                      #yes or no. An optional step for host-associated samples (eg. termite, human, plant...)
threads_bowtie2: 100                                                    #The number of threads to run this process. To be adjusted according to your hardware
host_genomes_directrory: /path/to/working/host_genomes/                 #the directory where the host genome is stored
host_genomes: /path/to/working/host_genomes/host_genomes.fa             #A fasta file containing the DNA sequences of the host genome(s)
threads_samtools: 100                                                   #The number of threads to run this process. To be adjusted according to your hardware
ressources_host_filtering: 500                                          #Memory according to tools need

################
### Assembly ###
################
threads_megahit: 150                                                   #The number of threads to run this process. To be adjusted according to your hardware
min_contig_len: 1000                                                   #Minimum length (in bp) of the assembled contigs
k_list: "21,31,41,51,61,71,81,91,99,109,119"                           #Kmer size (for further details, see the megahit documentation)
ressources_megahit: 500                                                #Memory according to tools need

###############
### Binning ###
###############
threads_bwa: 150                                                       #The number of threads to run this process. To be adjusted according to your hardware
ressources_bwa: 500                                                    #Memory according to tools need
threads_samtools: 150                                                  #The number of threads to run this process. To be adjusted according to your hardware
ressources_samtools: 500                                               #Memory according to tools need
seed: 19860615                                                         #Seed number for reproducible results
threads_metabat: 150                                                   #The number of threads to run this process. To be adjusted according to your hardware
minContig: 2500                                                        #Minimum length (in bp) of the contigs
ressources_binning: 500                                                #Memory according to tools need

####################
### Bins quality ###
####################
threads_checkm: 150                                                    #The number of threads to run this process. To be adjusted according to your hardware
ressources_checkm: 500                                                 #Memory according to tools need

######################
### Classification ###
######################
GTDB_data_ref: /path/to/downloaded/GTDB                                #Path to uncompressed GTDB-Tk reference data (GTDB)
threads_gtdb: 16                                                       #The number of threads to run this process. To be adjusted according to your hardware
ressources_gtdb: 500                                                   #Memory according to tools need

##################
### Abundances ###
##################
threads_coverM: 16                                                     #The number of threads to run this process. To be adjusted according to your hardware
ressources_coverM: 500                                                 #Memory according to tools need
```
# Run SnakeMAGs
If you are using a workstation with Ubuntu (tested on Ubuntu 22.04):
```{bash}
snakemake --cores 30 --snakefile SnakeMAGs.smk --use-conda --conda-prefix /path/to/SnakeMAGs_conda_env/ --configfile /path/to/config.yaml --keep-going --latency-wait 180
```

If you are working on a cluster with Slurm (tested with version 18.08.7):
```{bash}
snakemake --snakefile SnakeMAGs.smk --cluster 'sbatch -p <cluster_partition> --mem <memory> -c <cores> -o "cluster_logs/{wildcards}.{rule}.{jobid}.out" -e "cluster_logs/{wildcards}.{rule}.{jobid}.err" ' --jobs <nbr_of_parallel_jobs> --use-conda --conda-frontend conda --conda-prefix /path/to/SnakeMAGs_conda_env/ --jobname "{rule}.{wildcards}.{jobid}" --latency-wait 180 --configfile /path/to/config.yaml --keep-going
```

# Citations

If you use SnakeMAGs, please cite:
> Nachida Tadrent, Franck Dedeine, Vincent Hervé (In preparation). SnakeMAGs: a simple, efficient, flexible and scalable workflow to reconstruct prokaryotic genomes from metagenomes.

Please also cite the dependencies:
- [illumina-utils](https://doi.org/10.1371/journal.pone.0066643) : Murat Eren, A., Vineis, J. H., Morrison, H. G., & Sogin, M. L. (2013). A Filtering Method to Generate High Quality Short Reads Using Illumina Paired-End Technology. *PloS ONE*, 8(6), e66643.
- [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) : Bolger, A. M., Lohse, M., & Usadel, B. (2014). Genome analysis Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15), 2114-2120.
- [Bowtie2](https://doi.org/10.1038/nmeth.1923) : Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357–359.
- [SAMtools](https://doi.org/10.1093/bioinformatics/btp352) : Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. 
- [BEDtools](https://doi.org/10.1093/bioinformatics/btq033) : Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26(6), 841–842.
- [MEGAHIT](https://doi.org/10.1093/bioinformatics/btv033) : Li, D., Liu, C. M., Luo, R., Sadakane, K., & Lam, T. W. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31(10), 1674–1676. 
- [bwa](https://doi.org/10.1093/bioinformatics/btp324) : Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*, 25(14), 1754–1760. 
- [MetaBAT2](https://doi.org/10.7717/peerj.7359) : Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: An adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*, 2019(7), 1–13. 
- [CheckM](https://doi.org/10.1101/gr.186072.114) : Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*, 25(7), 1043–1055. 
- [GTDB-Tk](https://doi.org/10.1093/bioinformatics/btz848) : Chaumeil, P. A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2020). GTDB-Tk: A toolkit to classify genomes with the genome taxonomy database. *Bioinformatics*, 36(6), 1925–1927. 
- [CoverM](https://github.com/wwood/CoverM)

# License
This project is licensed under the CeCILL License - see the [LICENSE](https://github.com/Nachida08/SnakeMAGs/blob/main/LICENCE) file for details.

Developed by Nachida Tadrent at the Insect Biology Research Institute ([IRBI](https://irbi.univ-tours.fr/))
