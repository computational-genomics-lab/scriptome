import luigi
import os
import time
import subprocess
from tasks.rnaSeq.index_genome import indexGenome
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import filtlong
from tasks.readCleaning.reFormatReads import reformat
from tasks.annotation.make_tx_to_gene import makeTx2Gene

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

class indexTranscript(luigi.Task):

	project_name=luigi.Parameter(default="RNASeqAnalysis")
	genome_name = GlobalParameter().genome_name
	genome_suffix=GlobalParameter().genome_suffix
	read_library_type =GlobalParameter().read_library_type
	organism_domain = GlobalParameter().organism_domain
	adapter=GlobalParameter().adapter
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	quant_method = luigi.ChoiceParameter(choices=["salmon", "kallisto"], var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF", "NA"], var_type=str)
	Salmon_Index_Parameter = luigi.Parameter(default="--type quasi -k 31")
	Kallisto_Index_Parameter = luigi.Parameter(default="-k 31")
	sampleName = luigi.Parameter()


	def requires(self):

		if (self.pre_process_reads == "yes"):
			return [bbduk(sampleName=self.sampleName,read_library_type=self.read_library_type),
				makeTx2Gene(annotation_file_type=self.annotation_file_type)]
		if (self.pre_process_reads == "no"):
			return [reformat(sampleName=self.sampleName,read_library_type=self.read_library_type),
				   makeTx2Gene(annotation_file_type=self.annotation_file_type)]
		

	def output(self):

		TranscriptIndexFolder = os.path.join(os.getcwd(), self.project_name, "alignment_free_dea","transcript_index", self.genome_name + "_" + self.quant_method +"_index" + "/")

		if (self.quant_method == "salmon") and (self.read_library_type == "pe"):
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
					'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}

		if (self.quant_method == "salmon") and (self.read_library_type == "se"):
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
					'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}

		if (self.quant_method == "kallisto") and (self.read_library_type == "pe"):
			return {'out': luigi.LocalTarget(TranscriptIndexFolder + "/" + "kallisto.index")}

		if (self.quant_method == "kallisto") and (self.read_library_type == "se"):
			return {'out': luigi.LocalTarget(TranscriptIndexFolder + "/" + "kallisto.index")}

	def run(self):
		TranscriptIndexFolder = os.path.join(os.getcwd(), self.project_name, "alignment_free_dea","transcript_index", self.genome_name + "_" + self.quant_method  + "_index" + "/")
		transcriptFolder = os.path.join(os.getcwd(), self.project_name, "alignment_free_dea","transcript_index", self.genome_name + "_transcriptome" + "/")
		genomeFolder = os.path.join(os.getcwd(), GlobalParameter().genome_dir + "/")



		cmd_remove_genome_index = "rm {genomeFolder}{genome_name}.{genome_suffix}.fai".format(genomeFolder=genomeFolder,
																							  genome_name=self.genome_name, 
																							  genome_suffix=self.genome_suffix)

		print ("****** NOW RUNNING COMMAND ******: " + cmd_remove_genome_index)
		print (run_cmd(cmd_remove_genome_index))


		cmd_run_salmon_index_pe = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {transcriptFolder}{genome_name}.ffn " \
								  "-i {TranscriptIndexFolder}" \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					genome_name=self.genome_name,
					transcriptFolder=transcriptFolder,
					genomeFolder=genomeFolder,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_salmon_index_se = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {transcriptFolder}{genome_name}.ffn " \
								  "-i {TranscriptIndexFolder}" \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					genome_name=self.genome_name,
					genomeFolder=genomeFolder,
					transcriptFolder=transcriptFolder,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_kallisto_index_pe = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
									"cd {TranscriptIndexFolder}; " \
									"kallisto index " \
									"--index=kallisto.index {transcriptFolder}{genome_name}.ffn " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					genome_name=self.genome_name,
					genomeFolder=genomeFolder,
					transcriptFolder=transcriptFolder,
					Kallisto_Index_Parameter=self.Kallisto_Index_Parameter)

		cmd_run_kallisto_index_se = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
									"cd {TranscriptIndexFolder}; " \
									"kallisto index " \
									"--index=kallisto.index {transcriptFolder}{genome_name}.ffn " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					genome_name=self.genome_name,
					genomeFolder=genomeFolder,
					transcriptFolder=transcriptFolder,
					Kallisto_Index_Parameter=self.Kallisto_Index_Parameter)

		if (self.quant_method == "salmon") and (self.read_library_type == "pe"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_pe)
			print (run_cmd(cmd_run_salmon_index_pe))


		if (self.quant_method == "salmon") and (self.read_library_type == "se"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_se)
			print (run_cmd(cmd_run_salmon_index_se))



		if (self.quant_method == "kallisto") and (self.read_library_type == "pe"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_index_pe)
			print (run_cmd(cmd_run_kallisto_index_pe))


		if (self.quant_method == "kallisto") and (self.read_library_type == "se"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_index_se)
			print (run_cmd(cmd_run_kallisto_index_se))

