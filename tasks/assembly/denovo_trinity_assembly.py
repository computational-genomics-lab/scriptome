import luigi
import os
import time
import subprocess
import pandas as pd
from luigi import Parameter

from tasks.assembly.kmergenie import kmergenie_formater_bbduk
from tasks.assembly.kmergenie import kmergenie_formater_reformat
from tasks.assembly.kmergenie import optimal_kmer
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import filtlong
from tasks.readCleaning.reFormatReads import reformat


def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
						 shell=True,
						 universal_newlines=True,
						 stdout=subprocess.PIPE,
						 executable='/bin/bash')
	output = p.communicate()[0]
	return output

######################################
#Function to prepare Trinity paired-end input
def prepare_trinity_pe_bbduk(input_file):
	with open(input_file) as ifh:
		sample_name_list = ifh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'

		read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")

		left_read_name_list = [ read_folder + x + left_read_name_suffix for x in sample_name_list]
		right_read_name_list =[ read_folder + x + right_read_name_suffix for x in sample_name_list]

		left_reads = ','.join(left_read_name_list)
		right_reads = ','.join(right_read_name_list)

		Trinity_PE_Input = "--left " + left_reads + " --right " + right_reads
		return Trinity_PE_Input

#Function to prepare Trinity Single-end input
def prepare_trinity_se_bbduk(input_file):
	with open(input_file) as ifh:
		sample_name_list = ifh.read().splitlines()
		read_name_suffix = '.fastq'
		read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_SE_Reads" + "/")
		read_name_list = [ read_folder + x + read_name_suffix for x in sample_name_list]
		reads = ','.join(read_name_list)
		Trinity_SE_Input = "--single " + reads + " "
		return Trinity_SE_Input
######################################
def prepare_trinity_pe_reformat(input_file):
	with open(input_file) as ifh:
		sample_name_list = ifh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'

		read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")

		left_read_name_list = [ read_folder + x + left_read_name_suffix for x in sample_name_list]
		right_read_name_list =[ read_folder + x + right_read_name_suffix for x in sample_name_list]

		left_reads = ','.join(left_read_name_list)
		right_reads = ','.join(right_read_name_list)

		Trinity_PE_Input = "--left " + left_reads + " --right " + right_reads
		return Trinity_PE_Input

#Function to prepare Trinity Single-end input
def prepare_trinity_se_reformat(input_file):
	with open(input_file) as ifh:
		sample_name_list = ifh.read().splitlines()
		read_name_suffix = '.fastq'
		read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_SE_Reads" + "/")
		read_name_list = [ read_folder + x + read_name_suffix for x in sample_name_list]
		reads = ','.join(read_name_list)
		Trinity_SE_Input = "--single " + reads + " "
		return Trinity_SE_Input
############################################
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

class trinity(luigi.Task):
	project_name=luigi.Parameter(default="RNASeqAnalysis")
	adapter = GlobalParameter().adapter
	organism_domain = GlobalParameter().organism_domain
	threads = GlobalParameter().threads
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	read_library_type = GlobalParameter().read_library_type


	def requires(self):
		if self.pre_process_reads=="yes" and self.read_library_type=="pe":
			return [bbduk(read_library_type="pe",
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		if self.pre_process_reads=="yes" and self.read_library_type=="se":
			return [bbduk(read_library_type="se",
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]
		if self.pre_process_reads=="no" and self.read_library_type=="pe":
			return [reformat(read_library_type="pe",
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		if self.pre_process_reads=="no" and self.read_library_type=="se":
			return [reformat(read_library_type=self.read_library_type,
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]

	def output(self):
		assembled_transcript_folder = os.path.join(os.getcwd(), "RNASeqAnalysis", "denovo_assembly","trinity" +"_"+ self.read_library_type )

		return {'out': luigi.LocalTarget(assembled_transcript_folder + "/"+ "Trinity.fasta")}



	def run(self):
		assembled_transcript_folder = os.path.join(os.getcwd(), "RNASeqAnalysis", "denovo_assembly","trinity" +"_"+ self.read_library_type )
		if self.read_library_type=="pe":
			input_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")
		if self.read_library_type=="se":
			input_sample_list = os.path.join(os.getcwd(), "sample_list", "se_samples.lst")


		if self.pre_process_reads == "yes" and self.read_library_type == "pe":
			ReadFolder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
			cmd_trinity_read_input = prepare_trinity_pe_bbduk(input_sample_list)

		if self.pre_process_reads == "yes" and self.read_library_type == "se":
			ReadFolder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_SE_Reads" + "/")
			cmd_trinity_read_input = prepare_trinity_se_bbduk(input_sample_list)

		if self.pre_process_reads == "no" and self.read_library_type == "pe":
			ReadFolder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
			cmd_trinity_read_input = prepare_trinity_pe_reformat(input_sample_list)
		if self.pre_process_reads == "no" and self.read_library_type == "se":
			ReadFolder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_SE_Reads" + "/")
			cmd_trinity_read_input = prepare_trinity_se_reformat(input_sample_list)

		cmd_run_trinity = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
							 "cd {ReadFolder}; " \
							 "Trinity --no_normalize_reads --seqType fq " \
							 "--min_contig_length 500 " \
							 "--max_memory 10G " \
							 " {cmd_trinity_read_input} " \
							 "--output {assembled_transcript_folder} " \
							 "--CPU 4 " \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					ReadFolder=ReadFolder,
					cmd_trinity_read_input=cmd_trinity_read_input)

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity)
		print (run_cmd(cmd_run_trinity))








