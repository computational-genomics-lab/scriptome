import luigi
import os
import time
import subprocess
import pandas as pd
from tasks.rnaSeq.generate_transcript_count_file import alignmentFreeQuant

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
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	feature_type=luigi.Parameter()
	adapter=luigi.Parameter()
	result_tag=luigi.Parameter()	
	

class alignmentFreeDEA(luigi.Task):
	adapter = GlobalParameter().adapter
	organism_domain = GlobalParameter().organism_domain
	read_library_type = GlobalParameter().read_library_type
	result_tag = GlobalParameter().result_tag
	feature_type = GlobalParameter().feature_type
	genome_name = GlobalParameter().genome_name

	# Local parameters
	project_name = luigi.Parameter(default="RNASeqAnalysis")
	quant_method = luigi.ChoiceParameter(choices=["salmon", "kallisto"], var_type=str)
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	attribute_type = luigi.ChoiceParameter(default="gene_id", description='''Specify attribute type in GTF annotation. 
													 string(=[gene_id])''', choices=["gene_id", "transcript_id"], var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF"], var_type=str)
	strand_type = luigi.ChoiceParameter(default="0", choices=['0', '1', '2'],
									   description='''perform strand-specific read counting. int([=0]unstranded) 
										OR [=1] stranded] OR [=2] reversely-stranded. default[=0]''')

	dea_method = luigi.ChoiceParameter(choices=["deseq2", "edger"], var_type=str)
	report_name = luigi.Parameter(default="DEA_Report")
	factor_of_intrest = luigi.Parameter(default="conditions", description="factor of intrest column of the target file (string [=condititions]). ")
	reference_condition = luigi.Parameter(default="control", description="reference biological condition.  (string [=control]")
	# target_file = luigi.Parameter(description="path to the design/target file. (string [=target.tsv]")
	alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
	p_adjust_method = luigi.ChoiceParameter(default="BH", description="p-value adjustment method.", choices=["BH", "BY"], var_type=str)
	fit_type = luigi.ChoiceParameter(default="local", description="mean-variance relationship.", choices=["parametric", "local", "mean"],
									 var_type=str)
	size_factor = luigi.ChoiceParameter(default="median", description="method to estimate the size factors.", choices=["median", "shorth"],
										var_type=str)

#Local Parameter


	def requires(self):
		return [alignmentFreeQuant(quant_method=self.quant_method,
						annotation_file_type=self.annotation_file_type,
						pre_process_reads=self.pre_process_reads)
			   ]

	def output(self):
		resultFolder = os.path.join(os.getcwd(), self.project_name,
									"alignment_free_dea",
									"DEAnalysis",
									self.dea_method + "_" + self.quant_method + "_" + self.result_tag + "_" + self.read_library_type + "/")
		#STAR Aligner Output
		return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}		





	def run(self):
		target_file = os.path.join(os.getcwd(), "sample_list", "target.tsv")
		# rmd_DESeq2File = os.path.expanduser(os.path.join(('~'), 'scriptome', 'tasks', 'utility', "PlotDESEQ2.Rmd"))
		# rmd_edgeRFile = os.path.expanduser(os.path.join(('~'), 'scriptome', 'tasks', 'utility', "PlotEDGER.Rmd"))

		tx2genefolder = os.path.join(os.getcwd(), self.project_name,
									 "alignment_free_dea",
									 "transcript_index",
						     		  self.genome_name + "_transcriptome" + "/")

		resultFolder = os.path.join(os.getcwd(), self.project_name,
										"alignment_free_dea",
										"DEAnalysis",
										self.dea_method + "_" + self.quant_method + "_" + self.result_tag + "_" + self.read_library_type + "/")

		transcriptQuantFolder = os.path.join(os.getcwd(),
												 self.project_name, "alignment_free_dea",
												 "ReadQuant",
												 self.quant_method + "_quant_" + self.read_library_type + "/")


		cmd_run_salmon_DESeq2 = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
									"cd {resultFolder};" \
									"salmon_DESeq2.r " \
									"-t {target_file} " \
									"-q {transcriptQuantFolder} " \
									"-G {tx2genefolder} " \
									"-v {factor_of_intrest} " \
									"-c {reference_condition} " \
									"-f {fit_type} " \
									"-a {alpha} " \
									"-p {p_adjust_method} " \
									"-l {size_factor} " \
									"-T $(which PlotDESEQ2.Rmd)" \
				.format(resultFolder=resultFolder,
						target_file=target_file,
						transcriptQuantFolder=transcriptQuantFolder,
						tx2genefolder=tx2genefolder,
						factor_of_intrest=self.factor_of_intrest,
						reference_condition=self.reference_condition,
						fit_type=self.fit_type,
						alpha=self.alpha,
						p_adjust_method=self.p_adjust_method,
						size_factor=self.size_factor)

		cmd_run_kallisto_DESeq2 = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
									  "cd {resultFolder};" \
									  "kallisto_DESeq2.r " \
									  "-t {target_file} " \
									  "-q {transcriptQuantFolder} " \
									  "-G {tx2genefolder} " \
									  "-v {factor_of_intrest} " \
									  "-c {reference_condition} " \
									  "-f {fit_type} " \
									  "-a {alpha} " \
									  "-p {p_adjust_method} " \
									  "-l {size_factor} " \
									  "-T $(which PlotDESEQ2.Rmd)" \
				.format(resultFolder=resultFolder,
						target_file=target_file,
						transcriptQuantFolder=transcriptQuantFolder,
						tx2genefolder=tx2genefolder,
						factor_of_intrest=self.factor_of_intrest,
						reference_condition=self.reference_condition,
						fit_type=self.fit_type,
						alpha=self.alpha,
						p_adjust_method=self.p_adjust_method,
						size_factor=self.size_factor)

		cmd_run_salmon_edgeR = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
								   "cd {resultFolder};" \
								   "salmon_edgeR.r " \
								   "-t {target_file} " \
								   "-q {transcriptQuantFolder} " \
								   "-G {tx2genefolder} " \
								   "-v {factor_of_intrest} " \
								   "-c {reference_condition} " \
								   "-a {alpha} " \
								   "-p {p_adjust_method} " \
								   "-T $(which PlotEDGER.Rmd)" \
				.format(resultFolder=resultFolder,
						target_file=target_file,
						transcriptQuantFolder=transcriptQuantFolder,
						tx2genefolder=tx2genefolder,
						factor_of_intrest=self.factor_of_intrest,
						reference_condition=self.reference_condition,
						fit_type=self.fit_type,
						alpha=self.alpha,
						p_adjust_method=self.p_adjust_method,
						size_factor=self.size_factor)

		cmd_run_kallisto_edgeR = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
									 "cd {resultFolder};" \
									 "kallisto_edgeR.r " \
									 "-t {target_file} " \
									 "-q {transcriptQuantFolder} " \
									 "-G {tx2genefolder} " \
									 "-v {factor_of_intrest} " \
									 "-c {reference_condition} " \
									 "-a {alpha} " \
									 "-p {p_adjust_method} " \
									 "-T $(which PlotEDGER.Rmd)" \
				.format(resultFolder=resultFolder,
						target_file=target_file,
						transcriptQuantFolder=transcriptQuantFolder,
						tx2genefolder=tx2genefolder,
						factor_of_intrest=self.factor_of_intrest,
						reference_condition=self.reference_condition,
						fit_type=self.fit_type,
						alpha=self.alpha,
						p_adjust_method=self.p_adjust_method,
						size_factor=self.size_factor)

			##########Commands
			# DESEQ2
		if all([self.read_library_type == "pe", self.quant_method == "salmon", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_DESeq2)
			print(run_cmd(cmd_run_salmon_DESeq2))

		if all([self.read_library_type == "se", self.quant_method == "salmon", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_DESeq2)
			print(run_cmd(cmd_run_salmon_DESeq2))

		if all([self.read_library_type == "pe", self.quant_method == "kallisto", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_DESeq2)
			print(run_cmd(cmd_run_kallisto_DESeq2))

		if all([self.read_library_type == "se", self.quant_method == "kallisto", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_DESeq2)
			print(run_cmd(cmd_run_kallisto_DESeq2))

			# EDGER

		if all([self.read_library_type == "pe", self.quant_method == "salmon", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_edgeR)
			print(run_cmd(cmd_run_salmon_edgeR))

		if all([self.read_library_type == "se", self.quant_method == "salmon", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_edgeR)
			print(run_cmd(cmd_run_salmon_edgeR))

		if all([self.read_library_type == "pe", self.quant_method == "kallisto", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_edgeR)
			print(run_cmd(cmd_run_kallisto_edgeR))

		if all([self.read_library_type == "se", self.quant_method == "kallisto", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_edgeR)
			print(run_cmd(cmd_run_kallisto_edgeR))
