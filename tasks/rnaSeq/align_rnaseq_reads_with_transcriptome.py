import luigi
import os
from tasks.assembly.genome_guided_trinity_assembly import ggTrinity

class GlobalParameter(luigi.Config):
	genome_suffix=luigi.Parameter()
	read_library_type=luigi.Parameter()
	organism_domain=luigi.Parameter()
	genome_name=luigi.Parameter()
	genome_dir=luigi.Parameter()
	transcriptome_dir=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	feature_type=luigi.Parameter()
	adapter=luigi.Parameter()

class mapReadToTans(luigi.Task):
	project_name = luigi.Parameter(default="RNASeqAnalysis")
	adapter = GlobalParameter().adapter
	read_library_type = GlobalParameter().read_library_type
	organism_domain = GlobalParameter().organism_domain
	threads = GlobalParameter().threads
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	genome_name = GlobalParameter().genome_name
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread", "star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF"], var_type=str)
	sampleName = luigi.Parameter()

	def requires(self):

		if self.read_library_type == "pe":
				return [ggTrinity(pre_process_reads=self.pre_process_reads,
								  project_name=self.project_name,
								   annotation_file_type=self.annotation_file_type,
								  read_library_type=self.read_library_type,
								   rnaseq_aligner=self.rnaseq_aligner,
								   sampleName=i)
						for i in [line.strip()
								  for line in
								  open(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))]]

		if self.read_library_type == "se":
				return [ggTrinity(pre_process_reads=self.pre_process_reads,
								  project_name=self.project_name,
								   annotation_file_type=self.annotation_file_type,
								  read_library_type=self.read_library_type,
								   rnaseq_aligner=self.rnaseq_aligner,
								   sampleName=i)
						for i in [line.strip()
								  for line in
								  open(os.path.join(os.getcwd(), "sample_list", "se_samples.lst"))]]

	def output(self):

		gg_trinity_map_folder = os.path.join(os.getcwd(), self.project_name, "rnaseq_transcript_alignment",
											 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

		# STAR OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		# HISAT2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		# BOWTIE2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		# DART OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		# SEGEMEHL OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		# SUBREAD OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_trinity_map_folder + "/" + self.genome_name + ".bam")}

	def run(self):
		if self.pre_process_reads == "no":
			cleanedReadFolder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
		if self.pre_process_reads == "yes":
			cleanedReadFolder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")

		refGenomeFolder = os.path.join(GlobalParameter().genome_dir + "/")

		genomeIndexFolder = os.path.join(os.getcwd(), self.project_name, "genome_index",
										 self.genome_name + "_" + self.rnaseq_aligner + "_index" + "/")
		mapFolder = os.path.join(os.getcwd(), self.project_name, "rnaseq_genome_alignment",
								 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

		qualimapFolder = os.path.join(os.getcwd(), self.project_name, "QualiMAP",
									  self.rnaseq_aligner + "_qualimap_" + self.read_library_type, self.sampleName + "/")

		##########################################################################################################
		# 1 SEGEMEHL rnaseq_aligner                                                                                          #
		##########################################################################################################

		cmd_run_segemehl_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t {threads} " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName, genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t {threads} " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					genome_name=self.genome_name, genome_suffix=self.genome_suffix,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									   "cd {mapFolder}; " \
									   "segemehl.x -t {threads} " \
									   "-i {genomeIndexFolder}index.idx " \
									   "-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
									   "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									   "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
									   "-S " \
									   "|samtools view -bS - | samtools sort " \
									   "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName, genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									   "cd {mapFolder}; " \
									   "segemehl.x -t {threads} " \
									   "-i {genomeIndexFolder}index.idx " \
									   "-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
									   "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									   "-S " \
									   "|samtools view -bS - | samtools sort " \
									   "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName, genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		################################################################################################################
		# 1. HISAT2 rnaseq_aligner
		################################################################################################################
		cmd_run_hisat2_pe_prok = "[-d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder};hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								 "-p {threads} " \
								 "--max-intronlen 20 " \
								 "--no-spliced-alignment " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								 "| samtools view -bS - | samtools sort " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					sampleName=self.sampleName,
					genome_name=self.genome_name)

		cmd_run_hisat2_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								 "-p {threads} " \
								 "--max-intronlen 20 " \
								 "--no-spliced-alignment " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "| samtools view -bS - | samtools sort " \
								 "-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)

		cmd_run_hisat2_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								"-p {threads} " \
								"-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								"-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								"| samtools view -bS - | samtools sort " \
								"-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)

		cmd_run_hisat2_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								"-p {threads} " \
								"-U {cleanedReadFolder}{sampleName}.fastq " \
								"| samtools view -bS - | samtools sort -m 5G " \
								"-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)

		################################################################################################################
		# 2. STAR rnaseq_aligner
		################################################################################################################
		cmd_run_star_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							   "cd {mapFolder}; " \
							   "STAR --runMode alignReads " \
							   "--alignIntronMax 1 " \
							   "--outSAMtype BAM " \
							   "SortedByCoordinate " \
							   "--limitBAMsortRAM 16000000000 " \
							   "--alignSJDBoverhangMin 999 " \
							   "--runThreadN {threads} " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							   "{cleanedReadFolder}{sampleName}_R2.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_star_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							   "cd {mapFolder}; " \
							   "STAR --runMode alignReads " \
							   "--alignIntronMax 1 " \
							   "--outSAMtype BAM " \
							   "SortedByCoordinate " \
							   "--limitBAMsortRAM 16000000000 " \
							   "--alignSJDBoverhangMin 999 " \
							   "--runThreadN {threads} " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_star_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "STAR --runMode alignReads " \
							  "--outSAMstrandField intronMotif " \
							  "--outSAMtype BAM " \
							  "SortedByCoordinate " \
							  "--limitBAMsortRAM 16000000000 " \
							  "--runThreadN {threads} " \
							  "--genomeDir {genomeIndexFolder} " \
							  "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							  "{cleanedReadFolder}{sampleName}_R2.fastq " \
							  "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					threads=self.threads,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_star_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "STAR --runMode alignReads " \
							  "--outSAMstrandField intronMotif " \
							  "--outSAMtype BAM " \
							  "SortedByCoordinate " \
							  "--limitBAMsortRAM 16000000000 " \
							  "--runThreadN {threads} " \
							  "--genomeDir {genomeIndexFolder} " \
							  "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							  "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_star_bam_rename = "mv {mapFolder}{sampleName}_Aligned.sortedByCoord.out.bam " \
							  "{mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder, sampleName=self.sampleName)

		##########################################################################################################
		# 3 DART rnaseq_aligner                                                                                          #
		##########################################################################################################

		cmd_run_dart_map_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "dart -intron 1 " \
								   "-t {threads} " \
								   "-i {genomeIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
								   "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_dart_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "dart -intron 1 " \
								   "-t {threads} " \
								   "-i {genomeIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_dart_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								  "cd {mapFolder}; " \
								  "dart -t {threads} " \
								  "-i {genomeIndexFolder}index " \
								  "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
								  "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								  "-j {mapFolder}{sampleName}_junctions.tab " \
								  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_dart_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								  "cd {mapFolder}; " \
								  "dart -t {threads} " \
								  "-i {genomeIndexFolder}index " \
								  "-f {cleanedReadFolder}{sampleName}.fastq " \
								  "-j {mapFolder}{sampleName}_junctions.tab " \
								  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		# Sort Unsorted BAM file and rempove the unsorted bam file

		cmd_run_dart_sort_bam = "cd {mapFolder}; samtools sort -@8 {sampleName}_unsorted.bam -o {sampleName}.bam".format(
			mapFolder=mapFolder, sampleName=self.sampleName)

		cmd_run_dart_remove_unsorted_bam = "cd {mapFolder}; rm {sampleName}_unsorted.bam".format(
			mapFolder=mapFolder, sampleName=self.sampleName)

		################################################################################################################
		# 4 Bowtie2 rnaseq_aligner
		################################################################################################################

		cmd_run_bowtie2_map_pe = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "bowtie2 -x {genomeIndexFolder}index " \
								 "-p {threads} " \
								 "-a " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq |samtools view -bS - | samtools " \
								 "sort " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_bowtie2_map_se = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "bowtie2 -x {genomeIndexFolder}index " \
								 "-p {threads} " \
								 "-a " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "|samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		################################################################################################################
		# 5 SUBREAD rnaseq_aligner
		################################################################################################################

		cmd_run_subread_map_pe = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "subread-align -i {genomeIndexFolder}index " \
								 "-T {threads} " \
								 "-t 0 " \
								 "-r {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-R {cleanedReadFolder}{sampleName}_R2.fastq " \
								 "--SAMoutput | samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_subread_map_se = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "subread-align -i {genomeIndexFolder}index " \
								 "-T {threads} " \
								 "-t 0 " \
								 "-r {cleanedReadFolder}{sampleName}.fastq " \
								 "--SAMoutput | samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		#######################################################################################################################
		# Mapping Quality Assessment using Qualimap
		########################################################################################################################
		cmd_run_qualimap = "[ -d {qualimapFolder} ] || mkdir -p {qualimapFolder}; " \
						   "cd {qualimapFolder}; " \
						   "qualimap bamqc " \
						   "-bam {mapFolder}{sampleName}.bam " \
						   "--java-mem-size={maxMemory}G " \
						   "-outdir {qualimapFolder} " \
						   "-outfile {sampleName}_QualiMap " \
						   "-outformat PDF:HTML" \
			.format(qualimapFolder=qualimapFolder,
					maxMemory=self.maxMemory,
					mapFolder=mapFolder,
					sampleName=self.sampleName)

		########################################################################################################################
		# Call rnaseq_aligner commands
		########################################################################################################################
		# Run Bowtie2: Only for prokaryotes
		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_aligner == "bowtie2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_pe)
			print(run_cmd(cmd_run_bowtie2_map_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_aligner == "bowtie2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_se)
			print(run_cmd(cmd_run_bowtie2_map_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		# Run Subraed
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_pe)
			print(run_cmd(cmd_run_subread_map_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_pe)
			print(run_cmd(cmd_run_subread_map_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_se)
			print(run_cmd(cmd_run_subread_map_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_se)
			print(run_cmd(cmd_run_subread_map_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		# Run Segmehl
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_euk)
			print(run_cmd(cmd_run_segemehl_map_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_euk)
			print(run_cmd(cmd_run_segemehl_map_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_prok)
			print(run_cmd(cmd_run_segemehl_map_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_prok)
			print(run_cmd(cmd_run_segemehl_map_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		# Run DART

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_euk)
			print(run_cmd(cmd_run_dart_map_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_se_euk)
			print(run_cmd(cmd_run_dart_map_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_prok)
			print(run_cmd(cmd_run_dart_map_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_se_prok)
			print(run_cmd(cmd_run_dart_map_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		# Run HISAT2
		####
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_euk)
			print(run_cmd(cmd_run_hisat2_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_euk)
			print(run_cmd(cmd_run_hisat2_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_prok)
			print(run_cmd(cmd_run_hisat2_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_prok)
			print(run_cmd(cmd_run_hisat2_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		# Run STAR

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_prok)
			print(run_cmd(cmd_run_star_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_prok)
			print(run_cmd(cmd_run_star_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_euk)
			print(run_cmd(cmd_run_star_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_euk)
			print(run_cmd(cmd_run_star_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print(run_cmd(cmd_run_qualimap))


#####################################################################################################
class mapReadsToGenome(luigi.Task):
	# Global Parameters
	project_name = luigi.Parameter(default="RNASeqAnalysis")
	read_library_type = GlobalParameter().read_library_type
	adapter = luigi.Parameter(default="./tasks/utility/adapters.fasta.gz")
	genome_name = GlobalParameter().genome_name
	organism_domain = GlobalParameter().organism_domain
	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory

	# Local Parameters
	# sampleName = luigi.Parameter()
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread", "star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF"], var_type=str)

	#####################################################################################################3

	def requires(self):
		if self.read_library_type == "pe":
			return [alignReads(pre_process_reads=self.pre_process_reads,
							   annotation_file_type=self.annotation_file_type,
							   rnaseq_aligner=self.rnaseq_aligner,
							   sampleName=i)
					for i in [line.strip()
							  for line in
							  open(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))]]

		if self.read_library_type == "se":
			return [alignReads(pre_process_reads=self.pre_process_reads,
							   annotation_file_type=self.annotation_file_type,
							   rnaseq_aligner=self.rnaseq_aligner,
							   sampleName=i)
					for i in [line.strip()
							  for line in
							  open(os.path.join(os.getcwd(), "sample_list", "se_samples.lst"))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(), "task_logs", 'task.align.read.to.genome.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('Read Alignment finished at {t}'.format(t=timestamp))