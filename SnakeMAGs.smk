################################
#######    SnakeMAGs     #######
################################

##### WORKING DIRECTORY #####

workdir: config["working_dir"]


##### SOME VARIABLES ######

WDIR = config["working_dir"]
HOST_GENOME = config['host_genome']
FASTQ_DIR = config["raw_fastq"]
SFX_1 = config["suffix_1"]
SFX_2 = config["suffix_2"]


###### WILDCARDS #####

SMP, = glob_wildcards(FASTQ_DIR+"{smp}"+SFX_1)

##### RULES #####


### CONJUNCTION RULES ("all") ###

rule all:
        input:
                expand("{smp}/Classification/{smp}_classification.final", smp=SMP),
                expand("{smp}/MAGs_abundances/{smp}_coverage.tsv", smp=SMP)

### TASK RULES ###
READS1 = FASTQ_DIR+"{smp}"+SFX_1
READS2 = FASTQ_DIR+"{smp}"+SFX_2

rule quality_filtering:
	input:
		fwd=READS1, rev=READS2
	output:
		quality_passed_R1="{smp}/QC_fq/quality_filtering/{smp}-QUALITY_PASSED_R1.fastq", quality_passed_R2="{smp}/QC_fq/quality_filtering/{smp}-QUALITY_PASSED_R2.fastq"
	conda: config['conda_env']+"IU.yaml"
	benchmark:
		"{smp}/benchmarks/quality_filtering.benchmark.txt"
	threads: config['threads_filter']
	params: email=config['email'], wdir=config['working_dir']
	resources: mem=config['resources_filter']
	log: "{smp}/logs/quality_filtering.log"
	shell:
		"""
		(mkdir -p {wildcards.smp}/QC_fq/quality_filtering/
		echo -e sample"\t"r1"\t"r2"\n"{wildcards.smp}"\t"{input.fwd}"\t"{input.rev} > {wildcards.smp}/QC_fq/quality_filtering/{wildcards.smp}_samples.txt
		iu-gen-configs {wildcards.smp}/QC_fq/quality_filtering/{wildcards.smp}_samples.txt --e-mail {params.email} -o {wildcards.smp}/QC_fq/quality_filtering/
		iu-filter-quality-minoche {wildcards.smp}/QC_fq/quality_filtering/{wildcards.smp}.ini --ignore-deflines) 2> {log}
		"""



rule adapter_trimming:
        input:
                quality_passed_R1="{smp}/QC_fq/quality_filtering/{smp}-QUALITY_PASSED_R1.fastq", quality_passed_R2="{smp}/QC_fq/quality_filtering/{smp}-QUALITY_PASSED_R2.fastq"
        output:
                trimmed_1="{smp}/QC_fq/adapter_trimming/{smp}_1.trimmed.fastq", un_trimmed_1="{smp}/QC_fq/adapter_trimming/{smp}_1un.trimmed.fastq", trimmed_2="{smp}/QC_fq/adapter_trimming/{smp}_2.trimmed.fastq", un_trimmed_2="{smp}/QC_fq/adapter_trimming/{smp}_2un.trimmed.fastq"
        conda: config['conda_env']+"TRIMMOMATIC.yaml"
	benchmark:
                "{smp}/benchmarks/adapter_trimming.benchmark.txt"
        params: wdir=config['working_dir'], adapters=config['adapters'], trim_params=config['trim_params']
        resources: mem=config['resources_trim']
        threads: config['threads_trim']
	log: "{smp}/logs/adapter_trimming.log"
        shell:
                """
                (mkdir -p {wildcards.smp}/QC_fq/adapter_trimming
                trimmomatic PE -threads {threads} -summary {wildcards.smp}/QC_fq/adapter_trimming/{wildcards.smp}-STATS.txt -trimlog {wildcards.smp}/QC_fq/adapter_trimming/{wildcards.smp}-log.txt {input.quality_passed_R1} {input.quality_passed_R2} {output.trimmed_1} {output.un_trimmed_1} {output.trimmed_2} {output.un_trimmed_2} ILLUMINACLIP:{params.adapters}:{params.trim_params}) 2> {log}
                """

if HOST_GENOME == 'yes':

	rule host_genomes_indexing:
		input:
			host_genomes=config['host_genomes']
		output:
			config["host_genomes_directory"]+"host_genomes_indexing.final"
		conda: config['conda_env']+"BOWTIE2.yaml"
		benchmark:
			config["host_genomes_directory"]+"benchmarks/host_reads_mapping.benchmark.txt"
		params: wdir=config['working_dir'], host_genomes=config['host_genomes'], host_genomes_directory=config['host_genomes_directory']
		threads: config['threads_bowtie2']
		resources: mem=config['resources_host_filtering']
		log: config["host_genomes_directory"]+"logs/host_reads_mapping.log"
		shell:
			"""
			(bowtie2-build --threads {threads} {input} {params.host_genomes_directory}/host_DB
			touch {output}) 2> {log}
                        """

	rule host_reads_mapping:
        	input:
                	trimmed_1="{smp}/QC_fq/adapter_trimming/{smp}_1.trimmed.fastq", trimmed_2="{smp}/QC_fq/adapter_trimming/{smp}_2.trimmed.fastq", index_final=config["host_genomes_directory"]+"host_genomes_indexing.final"
        	output:
                	sam="{smp}/QC_fq/host_filtering/{smp}_all.sam"
        	conda:  config['conda_env']+"BOWTIE2.yaml"
		benchmark:
                	"{smp}/benchmarks/host_reads_mapping.benchmark.txt"
        	params: wdir=config['working_dir'], host_genomes=config['host_genomes'], host_genomes_directory=config['host_genomes_directory']
        	threads: config['threads_bowtie2']
        	resources: mem=config['resources_host_filtering']
		log: "{smp}/logs/host_reads_mapping.log"
		shell:
                	"""
			(mkdir -p {wildcards.smp}/QC_fq/host_filtering
			bowtie2 -x {params.host_genomes_directory}/host_DB -1 {input.trimmed_1} -2 {input.trimmed_2} --un-conc {wildcards.smp}/QC_fq/host_filtering/unmapped.fq --al-conc {wildcards.smp}/QC_fq/host_filtering/mapped.fq -S {output.sam} --threads {threads}) 2> {log}
			"""
	rule sam2bam:
		input:
			sam="{smp}/QC_fq/host_filtering/{smp}_all.sam"
		output:
			bam="{smp}/QC_fq/host_filtering/{smp}_bothEndsUnmapped_sorted.bam"
		conda: config['conda_env']+"SAMTOOLS.yaml"
		benchmark:
			"{smp}/benchmarks/sam2bam.benchmark.txt"
		threads: config['threads_samtools']
		resources: mem=config['resources_host_filtering']
		log: "{smp}/logs/sam2bam.log"
		shell:
			"""
			(samtools view -@ 100 -bS {input.sam} > {wildcards.smp}/QC_fq/host_filtering/all.bam
			samtools view -b -@ 100 -f 12 -F 256 {wildcards.smp}/QC_fq/host_filtering/all.bam > {wildcards.smp}/QC_fq/host_filtering/bothEndsUnmapped.bam
			samtools sort -n {wildcards.smp}/QC_fq/host_filtering/bothEndsUnmapped.bam --threads {threads}  > {output.bam}) 2> {log}
			"""
	rule host_reads_removal:
		input:
			bam="{smp}/QC_fq/host_filtering/{smp}_bothEndsUnmapped_sorted.bam"
		output:
			fq1="{smp}/QC_fq/host_filtering/{smp}_host_removed_R1.fastq", fq2="{smp}/QC_fq/host_filtering/{smp}_host_removed_R2.fastq"
		conda: config['conda_env']+"BEDTOOLS.yaml"
		benchmark:
                        "{smp}/benchmarks/host_reads_removing.benchmark.txt"
		resources: mem=config['resources_host_filtering']
		log: "{smp}/logs/host_reads_removing.log"
		shell:
			"""
			(bedtools bamtofastq -i {input.bam} -fq {output.fq1} -fq2 {output.fq2}) 2> {log}
			"""
	rule assembly_reads:
		input:
			fq1="{smp}/QC_fq/host_filtering/{smp}_host_removed_R1.fastq", fq2="{smp}/QC_fq/host_filtering/{smp}_host_removed_R2.fastq"
		output:
			"{smp}/Assembly/{smp}.contigs.fa"
		conda: config['conda_env']+"MEGAHIT.yaml"
		benchmark:
			"{smp}/benchmarks/assembly_reads.benchmark.txt"
		threads: config['threads_megahit']
		resources: mem=config['resources_megahit']
		params: min_contig_len=config['min_contig_len'], k_list=config['k_list'], wdir=config['working_dir']
		log: "{smp}/logs/assembly_reads.log"
		shell:
			"""
			(mkdir -p {wildcards.smp}/Assembly/
			megahit -1 {input.fq1} -2 {input.fq2} --memory 0.95  --num-cpu-threads {threads} --out-dir {wildcards.smp}/Assembly/megahit_out --out-prefix {wildcards.smp} --min-contig-len {params.min_contig_len} --k-list {params.k_list}
			ln -s {params.wdir}{wildcards.smp}/Assembly/megahit_out/{wildcards.smp}.contigs.fa {wildcards.smp}/Assembly/{wildcards.smp}.contigs.fa) 2> {log}
			"""
	rule depth_file_first_step:
        	input:
                	contigs="{smp}/Assembly/{smp}.contigs.fa", clean_R1="{smp}/QC_fq/host_filtering/{smp}_host_removed_R1.fastq", clean_R2="{smp}/QC_fq/host_filtering/{smp}_host_removed_R2.fastq"
        	output:
                	sam="{smp}/Binning/{smp}.contigs.sam"
        	conda: config['conda_env']+"BWA.yaml"
        	benchmark:
                	"{smp}/benchmarks/depth_file_first_step.benchmark.txt"
        	threads: config['threads_bwa']
        	resources: mem=config['resources_bwa']
		params: wdir=config['working_dir']
		log: "{smp}/logs/depth_file_first_step.log"
        	shell:
                	"""
                	(mkdir -p {wildcards.smp}/Binning
			ln -sf {params.wdir}/{input.contigs} {wildcards.smp}/Binning/{wildcards.smp}.contigs.fa
                	bwa index {wildcards.smp}/Binning/{wildcards.smp}.contigs.fa
                	bwa mem -t {threads} {wildcards.smp}/Binning/{wildcards.smp}.contigs.fa {input.clean_R1} {input.clean_R2} > {output.sam}) 2> {log}
                	"""

else:
	rule assembly:
		input:
			fq1="{smp}/QC_fq/adapter_trimming/{smp}_1.trimmed.fastq", fq2="{smp}/QC_fq/adapter_trimming/{smp}_2.trimmed.fastq"
		output:
			"{smp}/Assembly/{smp}.contigs.fa"
		conda: config['conda_env']+"MEGAHIT.yaml"
		benchmark:
			"{smp}/benchmarks/assembly_reads.benchmark.txt"
		threads: config['threads_megahit']
		resources: mem=config['resources_megahit']
		params: min_contig_len=config['min_contig_len'], k_list=config['k_list'], wdir=config['working_dir']
		log: "{smp}/logs/assembly_reads.log"
		shell:
			"""
			(mkdir -p {wildcards.smp}/Assembly/
                        megahit -1 {input.fq1} -2 {input.fq2} --memory 0.95  --num-cpu-threads {threads} --out-dir {wildcards.smp}/Assembly/megahit_out --out-prefix {wildcards.smp} --min-contig-len {params.min_contig_len} --k-list {params.k_list}
                        ln -s {params.wdir}{wildcards.smp}/Assembly/megahit_out/{wildcards.smp}.contigs.fa {wildcards.smp}/Assembly/{wildcards.smp}.contigs.fa) 2> {log}
			"""


	rule depth_file_first_step:
		input:
			contigs="{smp}/Assembly/{smp}.contigs.fa", clean_R1="{smp}/QC_fq/adapter_trimming/{smp}_1.trimmed.fastq", clean_R2="{smp}/QC_fq/adapter_trimming/{smp}_2.trimmed.fastq"
		output:
			sam="{smp}/Binning/{smp}.contigs.sam"
		conda: config['conda_env']+"BWA.yaml"
		benchmark:
			"{smp}/benchmarks/depth_file_first_step.benchmark.txt"
		threads: config['threads_bwa']
		resources: mem=config['resources_bwa']
		params: wdir=config['working_dir']
		log: "{smp}/logs/depth_file_first_step.log"
		shell:
			"""
			(mkdir -p {wildcards.smp}/Binning
			ln -sf {params.wdir}/{input.contigs} {wildcards.smp}/Binning/{wildcards.smp}.contigs.fa
			bwa index {wildcards.smp}/Binning/{wildcards.smp}.contigs.fa
			bwa mem -t {threads} {wildcards.smp}/Binning/{wildcards.smp}.contigs.fa {input.clean_R1} {input.clean_R2} > {output.sam}) 2> {log}
			"""
rule depth_file_second_step:
	input:
		sam="{smp}/Binning/{smp}.contigs.sam"
	output:
		bam="{smp}/Binning/{smp}.sorted.contigs.bam"
	conda: config['conda_env']+"SAMTOOLS.yaml"
	benchmark:
		"{smp}/benchmarks/depth_file_second_step.benchmark.txt"
	threads: config['threads_samtools']
	resources: mem=config['resources_samtools']
	log: "{smp}/logs/depth_file_second_step.log"
	shell:
		"""
		(samtools view -b -@ {threads} {input.sam} -o {wildcards.smp}/Binning/{wildcards.smp}.contigs.bam
		samtools sort -@ {threads} {wildcards.smp}/Binning/{wildcards.smp}.contigs.bam -o {output.bam}) 2> {log}
		"""

rule depth_file:
	input:
		bam="{smp}/Binning/{smp}.sorted.contigs.bam"
	output:
		depth="{smp}/Binning/{smp}.depth.txt"
	conda: config["conda_env"]+"METABAT2.yaml"
	benchmark:
		"{smp}/benchmarks/depth_file.benchmark.txt"
	resources: mem=config['resources_binning']
	log: "{smp}/logs/depth_file.log"
	shell:
		"""
		(jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}) 2> {log}
		"""

rule binning:
        input:
                depth="{smp}/Binning/{smp}.depth.txt", contigs="{smp}/Assembly/{smp}.contigs.fa"
        output:
                "{smp}/Binning/{smp}_binning.final"
        conda: config["conda_env"]+"METABAT2.yaml"
        benchmark:
                "{smp}/benchmarks/binning.benchmark.txt"
        threads: config['threads_metabat']
        resources: mem=config['resources_binning']
        params: seed=config["seed"], minContig=config["minContig"]
	log: "{smp}/logs/binning.log"
        shell:
                """
                (mkdir -p {wildcards.smp}/Binning/bins
                metabat2 -i {input.contigs}  -a {input.depth} -o {wildcards.smp}/Binning/bins/bin --minContig {params.minContig} --numThreads {threads} --seed {params.seed} -v
                touch {wildcards.smp}/Binning/{wildcards.smp}_binning.final) 2> {log}
                """

rule bins_quality_checkM:
	input:
		"{smp}/Binning/{smp}_binning.final"
	output:
		"{smp}/Bins_quality/{smp}_checkM.final"
	conda: config["conda_env"]+"CHECKM.yaml"
	benchmark:
                "{smp}/benchmarks/checkm.benchmark.txt"
	threads: config['threads_checkm']
	resources: mem=config['resources_checkm']
	log: "{smp}/logs/checkm.log"
	shell:
		"""
		(mkdir -p {wildcards.smp}/Bins_quality/
		checkm lineage_wf -f {wildcards.smp}/Bins_quality/CheckM.txt -t {threads} -x fa {wildcards.smp}/Binning/bins/ {wildcards.smp}/Bins_quality
		checkm qa {wildcards.smp}/Bins_quality/lineage.ms {wildcards.smp}/Bins_quality/ --threads {threads} --out_format 2 --tab_table --file {wildcards.smp}/Bins_quality/quality_summary
		touch {wildcards.smp}/Bins_quality/{wildcards.smp}_checkM.final) 2> {log}
		"""

rule bins_quality_filtering:
        input:
                "{smp}/Bins_quality/{smp}_checkM.final"
        output:
                "{smp}/Bins_quality/{smp}_bins_filtering.final"
        benchmark:
                "{smp}/benchmarks/bins_filtering.benchmark.txt"
        params: completion=config['completion'], contamination=config['contamination'], PARKS_SCORE=config['parks_quality_score']
        log: "{smp}/logs/bins_filtering.log"
        shell:
                """
		mkdir -p {wildcards.smp}/Bins_quality/MAGs_checkM
		if [ {params.PARKS_SCORE} == 'yes' ]
		then
			awk 'BEGIN {{FS=" "}} {{if (($13-(5*$14))>=50) {{system("cp {wildcards.smp}/Binning/bins/"$1".fa {wildcards.smp}/Bins_quality/MAGs_checkM/"$1".fa")}}}}' {wildcards.smp}/Bins_quality/CheckM.txt
			touch {wildcards.smp}/Bins_quality/{wildcards.smp}_bins_filtering.final
		else
			awk 'BEGIN {{FS=" "}} {{if ($13>={params.completion} && $14<={params.contamination}) {{system("cp {wildcards.smp}/Binning/bins/"$1".fa {wildcards.smp}/Bins_quality/MAGs_checkM/"$1".fa")}}}}' {wildcards.smp}/Bins_quality/CheckM.txt
			touch {wildcards.smp}/Bins_quality/{wildcards.smp}_bins_filtering.final
		fi
		"""


rule bins_quality_gunc:
	input:
		"{smp}/Bins_quality/{smp}_bins_filtering.final"
	output:
		"{smp}/Bins_quality/{smp}_gunc.final"
	conda: config["conda_env"]+"GUNC.yaml"
	benchmark:
		"{smp}/benchmarks/GUNC.benchmark.txt"
	params: GUNC_DB=config['GUNC_db'], GUNC=config['gunc']
	threads: config['threads_gunc']
	resources: mem=config['resources_gunc']
	log: "{smp}/logs/gunc.log"
	shell:
		"""
		if [ {params.GUNC} == 'yes' ]
		then
			(mkdir -p {wildcards.smp}/Bins_quality/gunc_output {wildcards.smp}/Bins_quality/MAGs
			gunc run --db_file {params.GUNC_DB} --input_dir {wildcards.smp}/Bins_quality/MAGs_checkM --threads {threads} --out_dir {wildcards.smp}/Bins_quality/gunc_output
			awk 'BEGIN {{FS=" "}} {{if ($13=="True") {{system("cp {wildcards.smp}/Bins_quality/MAGs_checkM/"$1".fa {wildcards.smp}/Bins_quality/MAGs/"$1".fa")}}}}' {wildcards.smp}/Bins_quality/gunc_output/GUNC.progenomes_2.1.maxCSS_level.tsv
			touch {wildcards.smp}/Bins_quality/{wildcards.smp}_gunc.final) 2> {log}
		else
			mv {wildcards.smp}/Bins_quality/MAGs_checkM {wildcards.smp}/Bins_quality/MAGs
			touch {wildcards.smp}/Bins_quality/{wildcards.smp}_gunc.final
		fi
		"""

rule classification:
	input:
		"{smp}/Bins_quality/{smp}_gunc.final"
	output:
		"{smp}/Classification/{smp}_classification.final"
	conda: config['conda_env']+"GTDBTK.yaml"
	benchmark:
		"{smp}/benchmarks/classification.benchmark.txt"
	params: GTDB=config['GTDB_data_ref']
	threads: config['threads_gtdb']
	resources: mem=config['resources_gtdb']
	log: "{smp}/logs/classification.log"
	shell:
		"""
		(mkdir -p {wildcards.smp}/Classification/
		export GTDBTK_DATA_PATH={params.GTDB}
		gtdbtk classify_wf --cpus {threads} --pplacer_cpus {threads} --extension fa --genome_dir {wildcards.smp}/Bins_quality/MAGs --out_dir {wildcards.smp}/Classification/
		touch {wildcards.smp}/Classification/{wildcards.smp}_classification.final) 2> {log}
		"""

if HOST_GENOME == 'yes':

	rule MAGs_abundances:
		input:
			final="{smp}/Bins_quality/{smp}_gunc.final", fq1="{smp}/QC_fq/host_filtering/{smp}_host_removed_R1.fastq", fq2="{smp}/QC_fq/host_filtering/{smp}_host_removed_R2.fastq"
		output:
			"{smp}/MAGs_abundances/{smp}_coverage.tsv"
		conda: config['conda_env']+"COVERM.yaml"
		benchmark:
			"{smp}/benchmarks/abundances.benchmark.txt"
		params: wdir=config['working_dir']
		threads: config['threads_coverM']
		resources: mem=config['resources_coverM']
		log: "{smp}/logs/abundances.log"
		shell:
			"""
			(mkdir -p {wildcards.smp}/MAGs_abundances/
			coverm genome --coupled {input.fq1} {input.fq2} --genome-fasta-files {wildcards.smp}/Bins_quality/MAGs/*.fa -o {output} --threads {threads}) 2> {log}
			"""

else:
	rule MAGs_abundances:
		input:
			final="{smp}/Bins_quality/{smp}_gunc.final", fq1="{smp}/QC_fq/adapter_trimming/{smp}_1.trimmed.fastq", fq2="{smp}/QC_fq/adapter_trimming/{smp}_2.trimmed.fastq"
		output:
			"{smp}/MAGs_abundances/{smp}_coverage.tsv"
		conda: config['conda_env']+"COVERM.yaml"
		benchmark:
			"{smp}/benchmarks/abundances.benchmark.txt"
		params: wdir=config['working_dir']
		threads: config['threads_coverM']
		resources: mem=config['resources_coverM']
		log: "{smp}/logs/abundances.log"
		shell:
			"""
			(mkdir -p {wildcards.smp}/MAGs_abundances/
			coverm genome --coupled {input.fq1} {input.fq2} --genome-fasta-files {wildcards.smp}/Bins_quality/MAGs/*.fa -o {output} --threads {threads}) 2> {log}
			"""
