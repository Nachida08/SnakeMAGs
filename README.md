[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)


# About SnakeMAGs
SnakeMAGs is a workflow to reconstruct prokaryotic genomes from metagenomes. The main purpose of SnakeMAGs is to process Illumina data from raw reads to metagenome-assembled genomes (MAGs).
SnakeMAGs is efficient, easy to handle and flexible to different projects. The workflow is CeCILL licensed, implemented in Snakemake (run on multiple cores) and available for Linux.
SnakeMAGs performed eight main steps:
- Quality filtering of the reads
- Adapter trimming
- Filtering of the host sequences (optional)
- Assembly
- Binning
- Evaluation of the quality of the bins
- Classification of the MAGs
- Estimation of the relative abundance of the MAGs


![scheme of workflow](SnakeMAGs_schema.jpg?raw=true)

# How to use SnakeMAGs
## Install conda
The easiest way to install and run SnakeMAGs is to use [conda](https://www.anaconda.com/products/distribution). These package managers will help you to easily install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Install and activate Snakemake environment
Note: The workflow was developed with Snakemake 7.0.0
```
conda activate

# First, set up your channel priorities
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Then, create a new environment for the Snakemake version you require
conda create -n snakemake_7.0.0 snakemake=7.0.0

# And activate it
conda activate snakemake_7.0.0
```

Alternatively, you can also install Snakemake via mamba:
```
# If you do not have mamba yet on your machine, you can install it with:
conda install -n base -c conda-forge mamba

# Then you can install Snakemake
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake

# And activate it
conda activate snakemake

```

## SnakeMAGs executable
The easiest way to procure SnakeMAGs and its related files is to clone the repository using git:
```
git clone https://github.com/Nachida08/SnakeMAGs.git
```
Alternatively, you can download the relevant files:
```
wget https://github.com/Nachida08/SnakeMAGs/blob/main/SnakeMAGs.smk https://github.com/Nachida08/SnakeMAGs/blob/main/config.yaml
```

## SnakeMAGs input files
- Illumina paired-end reads in FASTQ.
- Adapter sequence file ([adapter.fa](https://github.com/Nachida08/SnakeMAGs/blob/main/adapters.fa)).
- Host genome sequences in FASTA (if host_genome: "yes"), in case you work with host-associated metagenomes (e.g. human gut metagenome).

## Download Genome Taxonomy Database (GTDB)
GTDB-Tk requires ~66G+ of external data (GTDB) that need to be downloaded and unarchived. Because this database is voluminous, we let you decide where you want to store it.
SnakeMAGs do not download automatically GTDB, you have to do it:

```
#Download the latest release (tested with release207 and also more recently with release214)
#Note: SnakeMAGs uses GTDBtk v2.1.0 and therefore require release 207 as minimum version. See https://ecogenomics.github.io/GTDBTk/installing/index.html#installing for details.
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz
#Decompress
tar -xzvf *tar.gz
#This will create a folder called release207_v2
```
All you have to do now is to indicate the path to the database folder (in our example, the folder is called release207_v2 but we also tested it recently with release214) in the config file, Classification section.

## Download the GUNC database (required if gunc: "yes")
GUNC accepts either a progenomes or GTDB based reference database. Both can be downloaded using the ```gunc download_db``` command. For our study we used the default proGenome-derived GUNC database. It requires less resources with similar performance.

```
conda activate
# Install and activate GUNC environment
conda create --prefix /path/to/gunc_env
conda install -c bioconda metabat2 --prefix /path/to/gunc_env
source activate /path/to/gunc_env

#Download the proGenome-derived GUNC database (tested with gunc_db_progenomes2.1)
#Note: SnakeMAGs uses GUNC v1.0.5
gunc download_db -db progenomes /path/to/GUNC_DB
```
All you have to do now is to indicate the path to the GUNC database file in the config file,  Bins quality section.

## Edit config file
You need to edit the config.yaml file. In particular, you need to set the correct paths: for the working directory, to specify where are your fastq files, where you want to place the conda environments (that will be created using the provided .yaml files available in [SnakeMAGs_conda_env directory](https://github.com/Nachida08/SnakeMAGs/tree/main/SnakeMAGs_conda_env)), where are the adapters, where is GTDB and optionally where is the GUNC database and where is your host genome reference.

Lastly, you need to allocate the proper computational resources (threads, memory) for each of the main steps. These can be optimized according to your hardware.



Here is an example of a config file:

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
suffix_1: "_R1.fastq"                                                    #Main type of suffix for forward reads file (eg. _1.fastq or _R1.fastq or _r1.fastq or _1.fq or _R1.fq or _r1.fq )
suffix_2: "_R2.fastq"                                                    #Main type of suffix for reverse reads file (eg. _2.fastq or _R2.fastq or _r2.fastq or _2.fq or _R2.fq or _r2.fq )

###########################
### Conda environnemnts ###
###########################

conda_env: "/path/to/SnakeMAGs_conda_env/"                               #Path to the provided SnakeMAGs_conda_env directory which contains the yaml file for each conda environment

#########################
### Quality filtering ###
#########################
email: name.surname@your-univ.com                                        #Your e-mail address
threads_filter: 10                                                       #The number of threads to run this process. To be adjusted according to your hardware
resources_filter: 150                                                    #Memory according to tools need (in GB)

########################
### Adapter trimming ###
########################
adapters: /path/to/working/directory/adapters.fa                         #A fasta file contanning a set of various Illumina adaptors (this file is provided and is also available on github)
trim_params: "2:40:15"                                                   #For further details, see the Trimmomatic documentation
threads_trim: 10                                                         #The number of threads to run this process. To be adjusted according to your hardware
resources_trim: 150                                                      #Memory according to tools need (in GB)

######################
### Host filtering ###
######################
host_genome: "yes"                                                      #yes or no. An optional step for host-associated samples (eg. termite, human, plant...)
threads_bowtie2: 50                                                     #The number of threads to run this process. To be adjusted according to your hardware
host_genomes_directory: /path/to/working/host_genomes/                  #the directory where the host genome is stored
host_genomes: /path/to/working/host_genomes/host_genomes.fa             #A fasta file containing the DNA sequences of the host genome(s)
threads_samtools: 50                                                    #The number of threads to run this process. To be adjusted according to your hardware
resources_host_filtering: 150                                           #Memory according to tools need (in GB)

################
### Assembly ###
################
threads_megahit: 50                                                    #The number of threads to run this process. To be adjusted according to your hardware
min_contig_len: 1000                                                   #Minimum length (in bp) of the assembled contigs
k_list: "21,31,41,51,61,71,81,91,99,109,119"                           #Kmer size (for further details, see the megahit documentation)
resources_megahit: 250                                                 #Memory according to tools need (in GB)

###############
### Binning ###
###############
threads_bwa: 50                                                        #The number of threads to run this process. To be adjusted according to your hardware
resources_bwa: 150                                                     #Memory according to tools need (in GB)
threads_samtools: 50                                                   #The number of threads to run this process. To be adjusted according to your hardware
resources_samtools: 150                                                #Memory according to tools need (in GB)
seed: 19860615                                                         #Seed number for reproducible results
threads_metabat: 50                                                    #The number of threads to run this process. To be adjusted according to your hardware
minContig: 2500                                                        #Minimum length (in bp) of the contigs
resources_binning: 250                                                 #Memory according to tools need (in GB)

####################
### Bins quality ###
####################
#checkM
threads_checkm: 50                                                    #The number of threads to run this process. To be adjusted according to your hardware
resources_checkm: 250                                                 #Memory according to tools need (in GB)
#bins_quality_filtering
completion: 50                                                        #The minimum completion rate of bins
contamination: 10                                                     #The maximum contamination rate of bins
parks_quality_score: "yes"                                            #yes or no. If yes bins are filtered according to the Parks quality score (completion-5*contamination >= 50)
#GUNC
gunc: "yes"                                                           #yes or no. An optional step to detect and discard chimeric and contaminated genomes using the GUNC tool
threads_gunc: 50                                                      #The number of threads to run this process. To be adjusted according to your hardware
resources_gunc: 250                                                   #Memory according to tools need (in GB)
GUNC_db: /path/to/GUNC_DB/gunc_db_progenomes2.1.dmnd                  #Path to the downloaded GUNC database (see the readme file)

######################
### Classification ###
######################
GTDB_data_ref: /path/to/downloaded/GTDB                                #Path to uncompressed GTDB-Tk reference data (GTDB)
threads_gtdb: 10                                                       #The number of threads to run this process. To be adjusted according to your hardware
resources_gtdb: 250                                                    #Memory according to tools need (in GB)

##################
### Abundances ###
##################
threads_coverM: 10                                                     #The number of threads to run this process. To be adjusted according to your hardware
resources_coverM: 150                                                  #Memory according to tools need (in GB)
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

If you are working on a cluster with SGE (tested with version 8.1.9):
```{bash}
snakemake --snakefile SnakeMAGs.smk --cluster "qsub -cwd -V -q <short.q/long.q> -pe thread {threads} -e cluster_logs/{rule}.e{jobid} -o cluster_logs/{rule}.o{jobid}" --jobs <nbr_of_parallel_jobs> --use-conda --conda-frontend conda --conda-prefix /path/to/SnakeMAGs_conda_env/ --jobname "{rule}.{wildcards}.{jobid}" --latency-wait 180 --configfile /path/to/config.yaml --keep-going
```


# Test
We provide you a small data set in the [test](https://github.com/Nachida08/SnakeMAGs/tree/main/test) directory which will allow you to validate your instalation and take your first steps with SnakeMAGs. This data set is a subset from [ZymoBiomics Mock Community](https://www.zymoresearch.com/blogs/blog/zymobiomics-microbial-standards-optimize-your-microbiomics-workflow) (250K reads) used in this tutoriel [metagenomics_tutorial](https://github.com/pjtorres/metagenomics_tutorial).

1. Before getting started make sure you have cloned the SnakeMAGs repository or you have downloaded all the necessary files (SnakeMAGs.smk, config.yaml, chr19.fa.gz, insub732_2_R1.fastq.gz, insub732_2_R2.fastq.gz). See the [SnakeMAGs executable](#snakemags-executable) section.
2. Unzip the fastq files and the host sequences file.
```
gunzip fastqs/insub732_2_R1.fastq.gz fastqs/insub732_2_R2.fastq.gz host_genomes/chr19.fa.gz
```
3. For better organisation put all the read files in the same directory (eg. fastqs) and the host sequences file in a separate directory (eg. host_genomes)
4. Edit the config file (see [Edit config file](#edit-config-file) section)
5. Run the test (see [Run SnakeMAGs](#run-snakemags) section)

Note: the analysis of these files took 1159.32 secondes to complete on a Ubuntu 22.04 LTS with an Intel(R) Xeon(R) Silver 4210 CPU @ 2.20GHz x 40 processor, 96GB of RAM.

# Genome reference for host reads filtering
For host-associated samples, one can remove host sequences from the metagenomic reads by mapping these reads against a reference genome. In the case of termite gut metagenomes, we are providing [here](https://zenodo.org/record/6908287#.YuAdFXZBx8M) the relevant files (fasta and index files) from termite genomes.

Upon request, we can help you to generate these files for your own reference genome and make them available to the community.

NB. These steps of mapping generate voluminous files such as .bam and .sam. Depending on your disk space, you might want to delete these files after use.


# Use case
During the test phase of the development of SnakeMAGs, we used this workflow to process 10 publicly available termite gut metagenomes generated by Illumina sequencing, to ultimately reconstruct prokaryotic MAGs. These metagenomes were retrieved from the NCBI database using the following accession numbers: SRR10402454; SRR14739927; SRR8296321; SRR8296327; SRR8296329; SRR8296337; SRR8296343; DRR097505; SRR7466794; SRR7466795. They come from five different studies: Waidele et al, 2019; Tokuda et al, 2018; Romero Victorica et al, 2020; Moreira et al, 2021; and Calusinska et al, 2020.

## Download the Illumina pair-end reads
We use fasterq-dump tool to extract data in FASTQ-format from SRA-accessions. It is a commandline-tool which offers a faster solution for downloading those large files.

```
# Install and activate sra-tools environment
## Note: For this study we used sra-tools 2.11.0

conda activate
conda install -c bioconda sra-tools
conda activate sra-tools

# Download fastqs in a single directory
mkdir raw_fastq
cd raw_fastq
fasterq-dump <SRA-accession> --threads <threads_nbr> --skip-technical --split-3
```

## Download Genome reference for host reads filtering
```
mkdir host_genomes
cd host_genomes
wget https://zenodo.org/record/6908287/files/termite_genomes.fasta.gz
gunzip termite_genomes.fasta.gz
```

## Edit the config file
See [Edit config file](#edit-config-file) section.

## Run SnakeMAGs
```
conda activate snakemake_7.0.0
mkdir cluster_logs
snakemake --snakefile SnakeMAGs.smk --cluster 'sbatch -p <cluster_partition> --mem <memory> -c <cores> -o "cluster_logs/{wildcards}.{rule}.{jobid}.out" -e "cluster_logs/{wildcards}.{rule}.{jobid}.err" ' --jobs <nbr_of_parallel_jobs> --use-conda --conda-frontend conda --conda-prefix /path/to/SnakeMAGs_conda_env/ --jobname "{rule}.{wildcards}.{jobid}" --latency-wait 180 --configfile /path/to/config.yaml --keep-going
```

## Study results
The MAGs reconstructed from each metagenome and their taxonomic classification are available in this [repository](https://doi.org/10.5281/zenodo.7661004).

# Citations

If you use SnakeMAGs, please cite:
> Tadrent N, Dedeine F and Hervé V. SnakeMAGs: a simple, efficient, flexible and scalable workflow to reconstruct prokaryotic genomes from metagenomes [version 2; peer review: 2 approved]. F1000Research 2023, 11:1522 (https://doi.org/10.12688/f1000research.128091.2)


Please also cite the dependencies:
- [Snakemake](https://doi.org/10.12688/f1000research.29032.2) : Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021) Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. *F1000Research* 2021, 10:33.
- [illumina-utils](https://doi.org/10.1371/journal.pone.0066643) : Murat Eren, A., Vineis, J. H., Morrison, H. G., & Sogin, M. L. (2013). A Filtering Method to Generate High Quality Short Reads Using Illumina Paired-End Technology. *PloS ONE*, 8(6), e66643.
- [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) : Bolger, A. M., Lohse, M., & Usadel, B. (2014). Genome analysis Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15), 2114-2120.
- [Bowtie2](https://doi.org/10.1038/nmeth.1923) : Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357–359.
- [SAMtools](https://doi.org/10.1093/bioinformatics/btp352) : Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079.
- [BEDtools](https://doi.org/10.1093/bioinformatics/btq033) : Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26(6), 841–842.
- [MEGAHIT](https://doi.org/10.1093/bioinformatics/btv033) : Li, D., Liu, C. M., Luo, R., Sadakane, K., & Lam, T. W. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31(10), 1674–1676.
- [bwa](https://doi.org/10.1093/bioinformatics/btp324) : Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*, 25(14), 1754–1760.
- [MetaBAT2](https://doi.org/10.7717/peerj.7359) : Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: An adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*, 2019(7), 1–13.
- [CheckM](https://doi.org/10.1101/gr.186072.114) : Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*, 25(7), 1043–1055.
- [GTDB-Tk](https://doi.org/10.1093/BIOINFORMATICS/BTAC672) : Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., Parks, D. H. (2022). GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics*.
- [CoverM](https://github.com/wwood/CoverM)
- [Waidele et al, 2019](https://doi.org/10.1101/526038) : Waidele, L., Korb, J., Voolstra, C. R., Dedeine, F., & Staubach, F. (2019). Ecological specificity of the metagenome in a set of lower termite species supports contribution of the microbiome to adaptation of the host. *Animal Microbiome*, 1(1), 1–13.
- [Tokuda et al, 2018](https://doi.org/10.1073/pnas.1810550115) : Tokuda, G., Mikaelyan, A., Fukui, C., Matsuura, Y., Watanabe, H., Fujishima, M., & Brune, A. (2018). Fiber-associated spirochetes are major agents of hemicellulose degradation in the hindgut of wood-feeding higher termites. *Proceedings of the National Academy of Sciences of the United States of America*, 115(51), E11996–E12004.
- [Romero Victorica et al, 2020](https://doi.org/10.1038/s41598-020-60850-5) : Romero Victorica, M., Soria, M. A., Batista-García, R. A., Ceja-Navarro, J. A., Vikram, S., Ortiz, M., Ontañon, O., Ghio, S., Martínez-Ávila, L., Quintero García, O. J., Etcheverry, C., Campos, E., Cowan, D., Arneodo, J., & Talia, P. M. (2020). Neotropical termite microbiomes as sources of novel plant cell wall degrading enzymes. *Scientific Reports*, 10(1), 1–14.
- [Moreira et al, 2021](https://doi.org/10.3389/fevo.2021.632590) : Moreira, E. A., Persinoti, G. F., Menezes, L. R., Paixão, D. A. A., Alvarez, T. M., Cairo, J. P. L. F., Squina, F. M., Costa-Leonardo, A. M., Rodrigues, A., Sillam-Dussès, D., & Arab, A. (2021). Complementary contribution of Fungi and Bacteria to lignocellulose digestion in the food stored by a neotropical higher termite. *Frontiers in Ecology and Evolution*, 9(April), 1–12.
- [Calusinska et al, 2020](https://doi.org/10.1038/s42003-020-1004-3) : Calusinska, M., Marynowska, M., Bertucci, M., Untereiner, B., Klimek, D., Goux, X., Sillam-Dussès, D., Gawron, P., Halder, R., Wilmes, P., Ferrer, P., Gerin, P., Roisin, Y., & Delfosse, P. (2020). Integrative omics analysis of the termite gut system adaptation to Miscanthus diet identifies lignocellulose degradation enzymes. *Communications Biology*, 3(1), 1–12.
- [Orakov et al, 2021](https://doi.org/10.1186/s13059-021-02393-0) : Orakov, A., Fullam, A., Coelho, L. P., Khedkar, S., Szklarczyk, D., Mende, D. R., Schmidt, T. S. B., & Bork, P. (2021). GUNC: detection of chimerism and contamination in prokaryotic genomes. *Genome Biology*, 22(1).
- [Parks et al, 2015](https://doi.org/10.1101/gr.186072.114) : Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*, 25(7), 1043–1055.
# License
This project is licensed under the CeCILL License - see the [LICENSE](https://github.com/Nachida08/SnakeMAGs/blob/main/LICENCE) file for details.

Developed by Nachida Tadrent at the Insect Biology Research Institute ([IRBI](https://irbi.univ-tours.fr/)), under the supervision of Franck Dedeine and Vincent Hervé.
