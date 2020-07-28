class mapReadToTans(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string ["
											  "=MyProject])")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	genome_name = luigi.Parameter()
	sampleName = luigi.Parameter()
	read_library_type = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")
	organism_domain = luigi.Parameter(description="organism_domain of the organism . (string [=prokaryote] OR [=eukaryote])")
	rnaseq_aligner = luigi.Parameter(description="""name of the rnaseq_aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for organism_domain eukaryote """)
	
	def requires(self):

		return [genomeGuidedTransAssembly(projectName=self.projectName,
								read_library_type=self.read_library_type,
								organism_domain = self.organism_domain,
								rnaseq_aligner = self.rnaseq_aligner,
								genome_name=self.genome_name,
								adapter=self.adapter)]


	def output(self):
		gg_transcript_map_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "ReadQuant",
										   "GG_Transcript_bowtie2_map_" + self.read_library_type + "/")

		return {'out1': luigi.LocalTarget(gg_transcript_map_folder + "/" + self.sampleName + ".bam")}

	def run(self):

		gg_transcript_map_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "ReadQuant",
										   "GG_Transcript_bowtie2_map_" + self.read_library_type + "/")

		gg_transcript_index_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "ReadQuant",
										   "GG_Transcript_Bowtie2_Index_" + self.read_library_type + "/")

		clean_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "InputReads" +"/")




		################################################################################################################
		#4 Bowtie2 rnaseq_aligner
		################################################################################################################

		cmd_run_bowtie2_map_pe = "[ -d {gg_transcript_map_folder} ] || mkdir -p {gg_transcript_map_folder}; " \
								 "cd {gg_transcript_map_folder}; " \
								 "bowtie2 -x {gg_transcript_index_folder}{genome_name} " \
								 "-p 2 " \
								 "-a " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq |samtools view -bS - | samtools " \
								 "sort " \
								 "-o {gg_transcript_map_folder}{sampleName}.bam " \
			.format(gg_transcript_map_folder=gg_transcript_map_folder,
					sampleName=self.sampleName,
					genome_name=self.genome_name,
					cleanedReadFolder=clean_read_folder,
					gg_transcript_index_folder=gg_transcript_index_folder)


		cmd_run_bowtie2_map_se = "[ -d {gg_transcript_map_folder} ] || mkdir -p {gg_transcript_map_folder}; " \
								 "cd {gg_transcript_map_folder}; " \
								 "bowtie2 -x {gg_transcript_index_folder}{genome_name} " \
								 "-p 2 " \
								 "-a " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "|samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {gg_transcript_map_folder}{sampleName}.bam " \
			.format(gg_transcript_map_folder=gg_transcript_map_folder,
					sampleName=self.sampleName,
					genome_name=self.genome_name,
					cleanedReadFolder=clean_read_folder,
					gg_transcript_index_folder=gg_transcript_index_folder)


		if self.read_library_type == "paired":
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_pe)
			print (run_cmd(cmd_run_bowtie2_map_pe))

		if self.read_library_type == "single":
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_se)
			print (run_cmd(cmd_run_bowtie2_map_se))

####################################################################################################

class mapReadToGGTansript(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string ["
											  "=MyProject])")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	genome_name = luigi.Parameter()

	#sampleName = luigi.Parameter()

	read_library_type = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	organism_domain = luigi.Parameter(description="organism_domain of the organism . (string [=prokaryote] OR [=eukaryote])")

	rnaseq_aligner = luigi.Parameter(description="""name of the rnaseq_aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for organism_domain eukaryote """)

	
	
	def requires(self):

		return [mapReadToTans(projectName=self.projectName,
								read_library_type=self.read_library_type,
								organism_domain = self.organism_domain,
								genome_name=self.genome_name,
								rnaseq_aligner = self.rnaseq_aligner,
								adapter=self.adapter,
								sampleName=i)
				for i in [line.strip()
						  for line in
								open (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))]]


	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))

