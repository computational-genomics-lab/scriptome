import luigi
import os
import time
import subprocess
import pandas as pd
from tasks.assembly.transcriptome_assembly import DTA

def run_cmd(cmd):
    p = subprocess.Popen(cmd, bufsize=-1,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
    output = p.communicate()[0]
    return output

def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)

createFolder("task_logs")

def prepare_corset_input(input_file):
    df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
    df1.columns = ["SampleName", "Condition"]
    df2 = df1.set_index("SampleName", drop = False)
    df_groups = df2.groupby('Condition')

    condition_group = ",".join(df1['Condition'])
    sample_group = ",".join(df1['SampleName'])
    CorSet_Input = '-g {} -n {}'.format(condition_group, sample_group)
    return CorSet_Input


class GlobalParameter(luigi.Config):
	genome_suffix=luigi.Parameter()
	read_library_type=luigi.Parameter()
	organism_domain=luigi.Parameter()
	genome_name=luigi.Parameter()
	transcriptome_dir=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	adapter=luigi.Parameter()

class indexDAT(luigi.Task):

	project_name=luigi.Parameter(default="RNASeqAnalysis")
	genome_name = GlobalParameter().genome_name
	read_library_type =GlobalParameter().read_library_type
	organism_domain = GlobalParameter().organism_domain
	adapter=GlobalParameter().adapter
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	#annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF", "NA"], var_type=str)
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	Salmon_Index_Parameter = luigi.Parameter(default="--type quasi -k 31")	
	#sampleName = luigi.Parameter()


	def requires(self):
		return [DTA(project_name=self.project_name,
					pre_process_reads=self.pre_process_reads,
					rnaseq_assembler=self.rnaseq_assembler)]
		

	def output(self):

		TranscriptIndexFolder = os.path.join(os.getcwd(), self.project_name,
											 "denovo_assembly_based_dea",
											 "transcript_index",
											 self.genome_name +  "_" + self.rnaseq_assembler +"_salmon_index" + "/")

		if self.rnaseq_assembler=="rockhopper":
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
				'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}

		if self.rnaseq_assembler=="spades":
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
				'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}
		if self.rnaseq_assembler=="trinity":
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
				'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}


		
	def run(self):

		if self.rnaseq_assembler == "trinity":
			assembled_transcript = os.path.join(os.getcwd(), "RNASeqAnalysis", "denovo_assembly","trinity" +"_"+ self.read_library_type + "/" + "Trinity.fasta")
		if self.rnaseq_assembler == "rockhopper":
			assembled_transcript = os.path.join(os.getcwd(), "RNASeqAnalysis", "denovo_assembly","rockhopper_" + self.read_library_type + "/" + "transcripts.fna")
		if self.rnaseq_assembler == "spades":
			assembled_transcript = os.path.join(os.getcwd(), "RNASeqAnalysis", "denovo_assembly", "spades_" + self.read_library_type + "/" + "transcripts.fasta")


		TranscriptIndexFolder = os.path.join(os.getcwd(), self.project_name, "denovo_assembly_based_dea","transcript_index", self.genome_name + "_"
											 + self.rnaseq_assembler +"_salmon_index" + "/")

		cmd_run_salmon_index_pe = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {assembled_transcript} " \
								  "-i {TranscriptIndexFolder} " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					assembled_transcript=assembled_transcript,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_salmon_index_se = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {assembled_transcript} " \
								  "-i {TranscriptIndexFolder} " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					assembled_transcript=assembled_transcript,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)


		if self.read_library_type == "pe":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_pe)
			print (run_cmd(cmd_run_salmon_index_pe))


		if self.read_library_type == "se":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_se)
			print (run_cmd(cmd_run_salmon_index_se))



##Run Salmom Quantification
class denovoQuant(luigi.Task):
	sampleName = luigi.Parameter()
	project_name=luigi.Parameter(default="RNASeqAnalysis")
	organism_domain = GlobalParameter().organism_domain
	read_library_type = GlobalParameter().read_library_type
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	adapter=GlobalParameter().adapter
	genome_name=GlobalParameter().genome_name
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	threads= GlobalParameter().threads
	memory=GlobalParameter().maxMemory


	def requires(self):
		return [indexDAT(project_name=self.project_name,
									pre_process_reads=self.pre_process_reads,
							        rnaseq_assembler=self.rnaseq_assembler)]

	def output(self):
		salmon_quant_folder = os.path.join(os.getcwd(), self.project_name, "denovo_assembly_based_dea", "ReadQuant",
										   self.rnaseq_assembler + "_denovo_salmon_quant_" + self.read_library_type + "/")

		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

	def run(self):

		salmon_quant_folder = os.path.join(os.getcwd(), self.project_name, "denovo_assembly_based_dea","ReadQuant",
										   self.rnaseq_assembler + "_denovo_salmon_quant_" + self.read_library_type + "/")

		salmon_index_folder = os.path.join(os.getcwd(), self.project_name, "denovo_assembly_based_dea","transcript_index",
											 self.genome_name + "_" + self.rnaseq_assembler + "_salmon_index" + "/")

		salmon_map_folder = os.path.join(os.getcwd(),  self.project_name,  "denovo_assembly_based_dea", "transcriptome", self.rnaseq_assembler +
										 "_denovo_salmon_map_" + self.read_library_type + "/")


		if self.pre_process_reads == "yes" and self.read_library_type == "pe":
			clean_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")

		if self.pre_process_reads == "yes" and self.read_library_type == "se":
			clean_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_SE_Reads" + "/")

		if self.pre_process_reads == "no" and self.read_library_type == "pe":
			clean_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")

		if self.pre_process_reads == "no" and self.read_library_type == "se":
			clean_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_SE_Reads" + "/")




		cmd_run_salmon_quant_pe = "[ -d {salmon_quant_folder} ] || mkdir -p {salmon_quant_folder}; " \
								  "[ -d {salmon_map_folder} ] || mkdir -p {salmon_map_folder}; " \
								  "cd {salmon_map_folder}; " \
								  "salmon quant --no-version-check -p {threads} " \
								  "-i {salmon_index_folder} " \
								  "-l A " \
								  "-1 {clean_read_folder}{sampleName}_R1.fastq " \
								  "-2 {clean_read_folder}{sampleName}_R2.fastq " \
								  "--dumpEq " \
								  "--output {salmon_quant_folder}{sampleName} " \
								  "--validateMappings " \
								  "--writeMappings | samtools view -bS - | samtools sort -m {memory}G " \
								  "-o {salmon_map_folder}{sampleName}.bam" \
			.format(salmon_quant_folder=salmon_quant_folder,
					salmon_map_folder=salmon_map_folder,
					sampleName=self.sampleName,
					clean_read_folder=clean_read_folder,
					threads=self.threads,
					salmon_index_folder=salmon_index_folder,
					memory=self.memory)

		cmd_run_salmon_quant_se = "[ -d {salmon_quant_folder} ] || mkdir -p {salmon_quant_folder}; " \
								  "mkdir -p {salmon_map_folder}; " \
								  "cd {salmon_map_folder};" \
								  "salmon quant --no-version-check -p {threads} " \
								  "--seqBias " \
								  "--gcBias " \
								  "-i {salmon_index_folder} " \
								  "-l A " \
								  "-r {clean_read_folder}{sampleName}.fastq " \
								  "--dumpEq " \
								  "--output {salmon_quant_folder}{sampleName} " \
								  "--validateMappings " \
								  "--writeMappings | samtools view -bS - | samtools sort -m {memory}G " \
								  "-o {salmon_map_folder}{sampleName}.bam" \
			.format(salmon_quant_folder=salmon_quant_folder,
					salmon_map_folder=salmon_map_folder,
					sampleName=self.sampleName,
					clean_read_folder=clean_read_folder,
					salmon_index_folder=salmon_index_folder,memory=self.memory,
					threads=self.threads)

		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))



		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))



		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))


		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))



		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))



		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))


##########################################################################################################################

class quantifyDAT(luigi.Task):

	project_name=luigi.Parameter(default="RNASeqAnalysis")
	organism_domain = GlobalParameter().organism_domain
	read_library_type = GlobalParameter().read_library_type
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	adapter=GlobalParameter().adapter
	genome_name=GlobalParameter().genome_name
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	threads= GlobalParameter().threads



	def requires(self):
		if self.read_library_type =="pe":
			return [denovoQuant(rnaseq_assembler=self.rnaseq_assembler,pre_process_reads=self.pre_process_reads,
							project_name=self.project_name,
							sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]

		if self.read_library_type =="se":
			return [denovoQuant(rnaseq_assembler=self.rnaseq_assembler,
							project_name=self.project_name,pre_process_reads=self.pre_process_reads,
							sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))

#####################################################################################################################################################3
# Run Corsetanalysis
class clusterDAT(luigi.Task):
	project_name=luigi.Parameter(default="RNASeqAnalysis")
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	organism_domain = GlobalParameter().organism_domain
	read_library_type = GlobalParameter().read_library_type
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	adapter=GlobalParameter().adapter

	def requires(self):

		return [
				[quantifyDAT(project_name=self.project_name,
											 pre_process_reads=self.pre_process_reads,
											 rnaseq_assembler=self.rnaseq_assembler)]
				]


	def output(self):
		corset_folder = os.path.join(os.getcwd(), self.project_name, "ReadQuant",
									 self.rnaseq_assembler + "_DenovoSalmonQuant_" + self.read_library_type, "Corset/" + "/")

		return {'out1': luigi.LocalTarget(corset_folder + "clusters.txt"),
				'out2': luigi.LocalTarget(corset_folder + "counts.txt")
				}


	def run(self):
		corset_folder = os.path.join(os.getcwd(), self.project_name, "denovo_assembly_based_dea", "ReadQuant",
										    + "denovo_" + self.rnaseq_assembler +  "_corset_" +   self.read_library_type + "/")

		salmon_quant_folder = os.path.join(os.getcwd(), self.project_name, "denovo_assembly_based_dea", "ReadQuant",
										   self.rnaseq_assembler + "_denovo_salmon_quant_" + self.read_library_type + "/")

		assembled_transcript_folder = os.path.join(os.getcwd(), "RNASeqAnalysis", "denovo_assembly", self.rnaseq_assembler + "_" + self.read_library_type)

		super_transcript_folder = os.path.join(os.getcwd(), "RNASeqAnalysis", "denovo_assembly",  self.rnaseq_assembler + "_" + self.read_library_type, "SuperTranscript/")

		# Paired end reads

		# Command to generate Rockhopper Paired-end input
		input_group_file = os.path.join(os.getcwd(),"sample_list", "group.tsv")

		# Command to generate Rockhopper Paired-end input
		cmd_corset_read_input = prepare_corset_input(input_group_file)

		cmd_run_corset = "[ -d  {corset_folder} ] || mkdir -p {corset_folder}; cd {corset_folder}; " \
						 "corset -D 99999999999  {corset_read_input} " \
						 "-i salmon_eq_classes {salmon_quant_folder}*/aux_info/eq_classes.txt " \
			.format(corset_folder=corset_folder,
					corset_read_input=cmd_corset_read_input,
					salmon_quant_folder=salmon_quant_folder)

		cmd_run_supertrans = "python $(which Lace.py) " \
							 "--cores 1 " \
							 "{assembled_transcript_folder}transcripts.fna " \
							 "{corset_folder}clusters.txt " \
							 "--outputDir {super_transcript_folder} " \
			.format(super_transcript_folder=super_transcript_folder, corset_folder=corset_folder,
					assembled_transcript_folder=assembled_transcript_folder)

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset)
		print(run_cmd(cmd_run_corset))

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_supertrans)
		print(run_cmd(cmd_run_supertrans))