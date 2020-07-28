import luigi
import math
from Bio import SeqIO
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
from tasks.rnaSeq.index_genome import indexGenome
from tasks.rnaSeq.align_rnaseq_reads_with_genome import alignReadsToGenome
from tasks.rnaSeq.align_rnaseq_reads_with_genome import alignReadSetsToGenome



def genomeSAindexNbases(genome_fasta):
    with open(genome_fasta) as genome:
        totalsize=0
        for rec in SeqIO.parse(genome, 'fasta'):
            totalsize = totalsize + len(rec)
        log2 = math.log(totalsize, 2.0)
        index = (log2/2) - 1
        gsanb = int(index)
        return gsanb

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
		Trinity_SE_Input = "--single" + reads + " "
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
		Trinity_SE_Input = "--single" + reads + " "
		return Trinity_SE_Input
############################################
def prepare_corset_input(input_file):
    df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
    df1.columns = ["SampleName", "Condition"]
    df2 = df1.set_index("SampleName", drop = False)
    df_groups = df2.groupby('Condition')

    condition_group = ",".join(df1['Condition'])
    sample_group = ",".join(df1['SampleName'])
    CorSet_Input = '-g {} -n {}'.format(condition_group, sample_group)
    return CorSet_Input
##############################################
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
	result_tag=luigi.Parameter()


class ggTrinity(luigi.Task):

	project_name=luigi.Parameter(default="RNASeqAnalysis")
	adapter = GlobalParameter().adapter
	read_library_type = GlobalParameter().read_library_type
	organism_domain = GlobalParameter().organism_domain
	threads = GlobalParameter().threads
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	genome_name=GlobalParameter().genome_name
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star","hisat2","dart", "segemehl","bowtie2"],var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF"], var_type=str)

	maxMemory = luigi.Parameter(description="Maximum Memory in GB",default="20")
	maxIntron = luigi.Parameter(description="Maximum Intron Length",default="2000")
	minContigLength = luigi.Parameter(description="Minimum Contig Length",default="200")
	threads = luigi.Parameter(description="Number of threads to be used", default="20")


	def requires(self):

		if self.read_library_type == "pe":
			return [alignReadSetsToGenome(project_name=self.project_name,pre_process_reads=self.pre_process_reads,
									   annotation_file_type=self.annotation_file_type,only_align="false",
									   rnaseq_aligner=self.rnaseq_aligner)]

		if self.read_library_type == "se":
			return [alignReadSetsToGenome(project_name=self.project_name,pre_process_reads=self.pre_process_reads,
									   annotation_file_type=self.annotation_file_type,only_align="false",
									   rnaseq_aligner=self.rnaseq_aligner)]


#Trinity --genome_guided_bam hexcentricum.bam --genome_guided_max_intron 10000 --max_memory 20G --CPU 10



	def output(self):
		gg_assembled_transcript_folder = os.path.join(os.getcwd(), self.project_name,
													  "genome_guided_assembly",
													  "trinity_" + self.rnaseq_aligner +"_" + self.read_library_type + "/")

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		# HISAT2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		# BOWTIE2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		# DART OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		# SEGEMEHL OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

	def run(self):
		mapFolder = os.path.join(os.getcwd(), self.project_name,
											 "genome_guided_assembly_based_dea",
											 "rnaseq_genome_alignment",
											 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")




		gg_assembled_transcript_folder = os.path.join(os.getcwd(), self.project_name,
													  "genome_guided_assembly",
													  "trinity_" + self.rnaseq_aligner +"_" + self.read_library_type + "/")


		cmd_run_trinity_pe = "[ -d  {gg_assembled_transcript_folder} ] || mkdir -p {gg_assembled_transcript_folder}; " \
							 "Trinity " \
							 "--min_contig_length {minContigLength} " \
							 "--max_memory {maxMemory}G " \
							 "--genome_guided_bam {mapFolder}{genome_name}.bam " \
							 "--genome_guided_max_intron {maxIntron} " \
							 "--output {gg_assembled_transcript_folder} " \
							 "--CPU {threads} " \
							 "--include_supertranscripts " \
							 "--full_cleanup " \
			.format(gg_assembled_transcript_folder=gg_assembled_transcript_folder,
					mapFolder=mapFolder,
					minContigLength=self.minContigLength,
					maxMemory=self.maxMemory,
					genome_name=self.genome_name,
					maxIntron=self.maxIntron,
					threads=self.threads)

		cmd_run_trinity_se = "[ -d  {gg_assembled_transcript_folder} ] || mkdir -p {gg_assembled_transcript_folder}; " \
							 "Trinity " \
							 "--min_contig_length {minContigLength} " \
							 "--max_memory {maxMemory}G " \
							 "--genome_guided_bam {mapFolder}{genome_name}.bam " \
							 "--genome_guided_max_intron {maxIntron} " \
							 "--output {gg_assembled_transcript_folder} " \
							 "--CPU {threads} " \
							 "--include_supertranscripts " \
							 "--full_cleanup " \
			.format(gg_assembled_transcript_folder=gg_assembled_transcript_folder,
					mapFolder=mapFolder,
					minContigLength=self.minContigLength,
					maxMemory=self.maxMemory,
					genome_name=self.genome_name,
					maxIntron=self.maxIntron,
					threads=self.threads)

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		# HISAT2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		# BOWTIE2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie22", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		# DART OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		# SEGEMEHL OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		#SUBREAD
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print(run_cmd(cmd_run_trinity_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print(run_cmd(cmd_run_trinity_se))

class indexGGAT(luigi.Task):
	project_name = luigi.Parameter(default="RNASeqAnalysis")
	genome_name = GlobalParameter().genome_name
	genome_suffix= GlobalParameter().genome_suffix
	genome_dir=GlobalParameter().genome_dir
	organism_domain = GlobalParameter().organism_domain
	read_library_type= GlobalParameter().read_library_type
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star","hisat2","dart", "segemehl","bowtie2"],var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF","GTF"],var_type=str)
	pre_process_reads=luigi.ChoiceParameter(choices=["yes","no"],var_type=str)
	threads = GlobalParameter().threads

	def requires(self):

		if self.read_library_type == "pe":
				return [ggTrinity(pre_process_reads=self.pre_process_reads,
								  project_name=self.project_name,
								  annotation_file_type=self.annotation_file_type,
								  rnaseq_aligner=self.rnaseq_aligner)]

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
		transcriptIndexFolder = os.path.join(os.getcwd(), self.project_name, "genome_guided_assembly_based_dea", "transcript_index",
											 self.genome_name + "_" + self.rnaseq_aligner + "_index" +
											 "/")

		if all([self.read_library_type== "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(transcriptIndexFolder + "SAindex")}

		if all([self.read_library_type== "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(transcriptIndexFolder + "SAindex")}

		if all([self.read_library_type== "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(transcriptIndexFolder + "SAindex")}

		if all([self.read_library_type== "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(transcriptIndexFolder + "SAindex")}

		if all([self.read_library_type== "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index" + ".1.ht2")}
		if all([self.read_library_type== "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index" + ".1.ht2")}
		if all([self.read_library_type== "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index" + ".1.ht2")}
		if all([self.read_library_type== "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index" + ".1.ht2")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "segemehl")]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index.idx")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "subread")]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index.files")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "dart")]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index.bwt")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "bowtie2")]):
			return {'out1': luigi.LocalTarget(transcriptIndexFolder + "index.1.bt2")}

	def run(self):
		refGenomeFolder = os.path.join(GlobalParameter().genome_dir + "/")

		transcriptomeFolder=os.path.join(os.getcwd(), self.project_name, "genome_guided_assembly",
													  "trinity_" + self.rnaseq_aligner + "_" + self.read_library_type + "/")

		transcriptIndexFolder = os.path.join(os.getcwd(), self.project_name, "genome_guided_assembly_based_dea","transcript_index",
											 self.genome_name + "_" + self.rnaseq_aligner + "_index" +
											 "/")

		transcriptFastaFile = "{transcriptomeFolder}Trinity-GG.fasta".format(transcriptomeFolder=transcriptomeFolder)

		cmd_gff2gtf = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder}; " \
					  "gffread -E  {refGenomeFolder}{genome_name}.gff -T " \
					  "-o {transcriptIndexFolder}{genome_name}.gtf " \
			.format(transcriptIndexFolder=transcriptIndexFolder,
					genome_name=self.genome_name,
					refGenomeFolder=refGenomeFolder)

		cmd_copy_gtf = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder}; " \
					   "cp {refGenomeFolder}{genome_name}.gtf {transcriptIndexFolder}{genome_name}.gtf " \
			.format(transcriptIndexFolder=transcriptIndexFolder,
					genome_name=self.genome_name,
					refGenomeFolder=refGenomeFolder)

		gsan = genomeSAindexNbases(transcriptFastaFile)

		cmd_run_star_index = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder}; STAR " \
							 "--runMode genomeGenerate " \
							 "--genomeSAindexNbases {gsan} " \
							 "--genomeFastaFiles {transcriptFastaFile} " \
							 "--genomeDir {transcriptIndexFolder} " \
			.format(transcriptIndexFolder=transcriptIndexFolder,
					transcriptFastaFile=transcriptFastaFile,
					gsan=gsan,
					refGenomeFolder=refGenomeFolder)

		cmd_run_bowtie2_index = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder};" \
								"bowtie2-build -f {transcriptFastaFile} {transcriptIndexFolder}index " \
			.format(transcriptIndexFolder=transcriptIndexFolder,
					transcriptFastaFile=transcriptFastaFile)

		cmd_run_dart_index = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder};" \
							 "bwt_index {transcriptFastaFile} {transcriptIndexFolder}index " \
			.format(transcriptIndexFolder=transcriptIndexFolder,
					transcriptFastaFile=transcriptFastaFile)

		cmd_run_hisat2_index = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder};" \
							   "hisat2-build {transcriptFastaFile} {transcriptIndexFolder}index " \
			.format(transcriptIndexFolder=transcriptIndexFolder,
					transcriptFastaFile=transcriptFastaFile)

		cmd_run_segemehl_index = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder};" \
								 "segemehl.x -x {transcriptIndexFolder}index.idx -d {transcriptFastaFile}  " \
			.format(transcriptIndexFolder=transcriptIndexFolder,transcriptFastaFile=transcriptFastaFile)

		cmd_run_subread_index = "[ -d  {transcriptIndexFolder} ] || mkdir -p {transcriptIndexFolder};" \
								"subread-buildindex -o {transcriptIndexFolder}index {transcriptFastaFile} " \
			.format(transcriptIndexFolder=transcriptIndexFolder,transcriptFastaFile=transcriptFastaFile)

		if (self.rnaseq_aligner == "star") and (self.annotation_file_type == "GTF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_index)
			print(run_cmd(cmd_run_star_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_copy_gtf)
			print(run_cmd(cmd_copy_gtf))
		if (self.rnaseq_aligner == "star") and (self.annotation_file_type == "GFF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_index)
			print(run_cmd(cmd_run_star_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "dart") and (self.annotation_file_type == "GTF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_index)
			print(run_cmd(cmd_run_dart_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_copy_gtf)
			print(run_cmd(cmd_copy_gtf))
		if (self.rnaseq_aligner == "dart") and (self.annotation_file_type == "GFF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_index)
			print(run_cmd(cmd_run_dart_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "bowtie2") and (self.annotation_file_type == "GTF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_index)
			print(run_cmd(cmd_run_bowtie2_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_copy_gtf)
			print(run_cmd(cmd_copy_gtf))
		if (self.rnaseq_aligner == "bowtie2") and (self.annotation_file_type == "GFF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_index)
			print(run_cmd(cmd_run_bowtie2_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "hisat2") and (self.annotation_file_type == "GTF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_index)
			print(run_cmd(cmd_run_hisat2_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_copy_gtf)
			print(run_cmd(cmd_copy_gtf))

		if (self.rnaseq_aligner == "hisat2") and (self.annotation_file_type == "GFF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_index)
			print(run_cmd(cmd_run_hisat2_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "segemehl") and (self.annotation_file_type == "GTF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_index)
			print(run_cmd(cmd_run_segemehl_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_copy_gtf)
			print(run_cmd(cmd_copy_gtf))
		if (self.rnaseq_aligner == "segemehl") and (self.annotation_file_type == "GFF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_index)
			print(run_cmd(cmd_run_segemehl_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "subread") and (self.annotation_file_type == "GTF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_index)
			print(run_cmd(cmd_run_subread_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_copy_gtf)
			print(run_cmd(cmd_copy_gtf))
		if (self.rnaseq_aligner == "subread") and (self.annotation_file_type == "GFF"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_index)
			print(run_cmd(cmd_run_subread_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))


class alingnReadToGGAT(luigi.Task):
	project_name = luigi.Parameter(default="RNASeqAnalysis")
	genome_name = GlobalParameter().genome_name
	genome_suffix = GlobalParameter().genome_suffix
	genome_dir = GlobalParameter().genome_dir
	organism_domain = GlobalParameter().organism_domain
	read_library_type = GlobalParameter().read_library_type
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread", "star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF"], var_type=str)
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	sampleName=luigi.Parameter()
	threads = GlobalParameter().threads

	def requires(self):

		return [indexGGAT(project_name=self.project_name,
									   rnaseq_aligner=self.rnaseq_aligner,
									   pre_process_reads=self.pre_process_reads,
									   annotation_file_type=self.annotation_file_type)]

	def output(self):
		mapFolder = os.path.join(os.getcwd(), self.project_name,
								 "genome_guided_assembly_based_dea",
								 "rnaseq_transcript_alignment",
								 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

		# STAR OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		# HISAT2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		# BOWTIE2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie22", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		# DART OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		# SEGEMEHL OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		# SUBREAD OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}


	def run(self):
		transcriptomeFolder = os.path.join(os.getcwd(), self.project_name, "genome_guided_assembly",
										   "trinity_" + self.rnaseq_aligner + "_" + self.read_library_type + "/")

		transcriptFastaFile = "{transcriptomeFolder}Trinity-GG.fasta".format(transcriptomeFolder=transcriptomeFolder)


		mapFolder = os.path.join(os.getcwd(), self.project_name,
								 "genome_guided_assembly_based_dea",
								 "rnaseq_transcript_alignment",
								 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

		transcriptIndexFolder = os.path.join(os.getcwd(), self.project_name,
										 "genome_guided_assembly_based_dea",
										 "transcript_index",
										 self.genome_name + "_" + self.rnaseq_aligner + "_index" +	 "/")



		if self.pre_process_reads == "yes" and self.read_library_type == "pe":
			cleanedReadFolder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
		if self.pre_process_reads == "yes" and self.read_library_type == "se":
			cleanedReadFolder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_SE_Reads" + "/")
		if self.pre_process_reads == "no" and self.read_library_type == "pe":
			cleanedReadFolder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
		if self.pre_process_reads == "no" and self.read_library_type == "se":
			cleanedReadFolder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_SE_Reads" + "/")

		refGenomeFolder = os.path.join(GlobalParameter().genome_dir + "/")



		##########################################################################################################
		# 1 SEGEMEHL rnaseq_aligner                                                                                          #
		##########################################################################################################

		cmd_run_segemehl_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t {threads} " \
									  "-i {transcriptIndexFolder}index.idx " \
									  "-d {transcriptFastaFile} " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads = GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder,
					transcriptFastaFile=transcriptFastaFile)

		cmd_run_segemehl_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t {threads} " \
									  "-i {transcriptIndexFolder}index.idx " \
									  "-d {transcriptFastaFile} " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads = GlobalParameter().threads,
					sampleName=self.sampleName,
					transcriptFastaFile=transcriptFastaFile,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_segemehl_map_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									   "cd {mapFolder}; " \
									   "segemehl.x -t {threads} " \
									   "-i {transcriptIndexFolder}index.idx " \
									   "-d {transcriptFastaFile} " \
									   "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									   "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
									   "-S " \
									   "|samtools view -bS - | samtools sort " \
									   "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads = GlobalParameter().threads,
					sampleName=self.sampleName,
					transcriptFastaFile=transcriptFastaFile,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_segemehl_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									   "cd {mapFolder}; " \
									   "segemehl.x -t {threads} " \
									   "-i {transcriptIndexFolder}index.idx " \
									   "-d {transcriptFastaFile} " \
									   "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									   "-S " \
									   "|samtools view -bS - | samtools sort " \
									   "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder,
					transcriptFastaFile=transcriptFastaFile)

		################################################################################################################
		# 1. HISAT2 rnaseq_aligner
		################################################################################################################
		cmd_run_hisat2_pe_prok = "[-d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder};hisat2 --dta -x {transcriptIndexFolder}index " \
								 "-p {threads} " \
								 "--max-intronlen 20 " \
								 "--no-spliced-alignment " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								 "| samtools view -bS - | samtools sort " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads = GlobalParameter().threads,
					transcriptIndexFolder=transcriptIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					sampleName=self.sampleName,
					genome_name=self.genome_name)

		cmd_run_hisat2_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; hisat2 --dta -x {transcriptIndexFolder}index " \
								 "-p {threads} " \
								 "--max-intronlen 20 " \
								 "--no-spliced-alignment " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "| samtools view -bS - | samtools sort " \
								 "-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads = GlobalParameter().threads,
					transcriptIndexFolder=transcriptIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)

		cmd_run_hisat2_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder}; hisat2 --dta -x {transcriptIndexFolder}index " \
								"-p {threads} " \
								"-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								"-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								"| samtools view -bS - | samtools sort " \
								"-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					transcriptIndexFolder=transcriptIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)

		cmd_run_hisat2_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder}; hisat2 --dta -x {transcriptIndexFolder}index " \
								"-p {threads} " \
								"-U {cleanedReadFolder}{sampleName}.fastq " \
								"| samtools view -bS - | samtools sort -m 5G " \
								"-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					transcriptIndexFolder=transcriptIndexFolder,
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
							   "--genomeDir {transcriptIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							   "{cleanedReadFolder}{sampleName}_R2.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_star_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							   "cd {mapFolder}; " \
							   "STAR --runMode alignReads " \
							   "--alignIntronMax 1 " \
							   "--outSAMtype BAM " \
							   "SortedByCoordinate " \
							   "--limitBAMsortRAM 16000000000 " \
							   "--alignSJDBoverhangMin 999 " \
							   "--runThreadN {threads} " \
							   "--genomeDir {transcriptIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_star_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "STAR --runMode alignReads " \
							  "--outSAMstrandField intronMotif " \
							  "--outSAMtype BAM " \
							  "SortedByCoordinate " \
							  "--limitBAMsortRAM 16000000000 " \
							  "--runThreadN {threads} " \
							  "--genomeDir {transcriptIndexFolder} " \
							  "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							  "{cleanedReadFolder}{sampleName}_R2.fastq " \
							  "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					threads=GlobalParameter().threads,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_star_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "STAR --runMode alignReads " \
							  "--outSAMstrandField intronMotif " \
							  "--outSAMtype BAM " \
							  "SortedByCoordinate " \
							  "--limitBAMsortRAM 16000000000 " \
							  "--runThreadN {threads} " \
							  "--genomeDir {transcriptIndexFolder} " \
							  "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							  "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

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
								   "-i {transcriptIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
								   "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_dart_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "dart -intron 1 " \
								   "-t {threads} " \
								   "-i {transcriptIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_dart_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								  "cd {mapFolder}; " \
								  "dart -t {threads} " \
								  "-i {transcriptIndexFolder}index " \
								  "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
								  "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								  "-j {mapFolder}{sampleName}_junctions.tab " \
								  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_dart_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								  "cd {mapFolder}; " \
								  "dart -t {threads} " \
								  "-i {transcriptIndexFolder}index " \
								  "-f {cleanedReadFolder}{sampleName}.fastq " \
								  "-j {mapFolder}{sampleName}_junctions.tab " \
								  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

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
								 "bowtie2 -x {transcriptIndexFolder}index " \
								 "-p {threads} " \
								 "-a " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq |samtools view -bS - | samtools " \
								 "sort " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_bowtie2_map_se = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "bowtie2 -x {transcriptIndexFolder}index " \
								 "-p {threads} " \
								 "-a " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "|samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		################################################################################################################
		# 5 SUBREAD rnaseq_aligner
		################################################################################################################

		cmd_run_subread_map_pe = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "subread-align -i {transcriptIndexFolder}index " \
								 "-T {threads} " \
								 "-t 0 " \
								 "-r {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-R {cleanedReadFolder}{sampleName}_R2.fastq " \
								 "--SAMoutput | samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		cmd_run_subread_map_se = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "subread-align -i {transcriptIndexFolder}index " \
								 "-T {threads} " \
								 "-t 0 " \
								 "-r {cleanedReadFolder}{sampleName}.fastq " \
								 "--SAMoutput | samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					transcriptIndexFolder=transcriptIndexFolder)

		#######################################################################################################################


		########################################################################################################################
		# Call rnaseq_aligner commands
		########################################################################################################################
		# Run Bowtie2: Only for prokaryotes
		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_aligner == "bowtie2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_pe)
			print(run_cmd(cmd_run_bowtie2_map_pe))



		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_aligner == "bowtie2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_se)
			print(run_cmd(cmd_run_bowtie2_map_se))



		# Run Subraed
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_pe)
			print(run_cmd(cmd_run_subread_map_pe))


		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_pe)
			print(run_cmd(cmd_run_subread_map_pe))



		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_se)
			print(run_cmd(cmd_run_subread_map_se))


		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_se)
			print(run_cmd(cmd_run_subread_map_se))



		# Run Segmehl
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_euk)
			print(run_cmd(cmd_run_segemehl_map_pe_euk))



		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_euk)
			print(run_cmd(cmd_run_segemehl_map_se_euk))


		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_prok)
			print(run_cmd(cmd_run_segemehl_map_pe_prok))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_prok)
			print(run_cmd(cmd_run_segemehl_map_se_prok))
	# Run DART

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_euk)
			print(run_cmd(cmd_run_dart_map_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))


			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_se_euk)
			print(run_cmd(cmd_run_dart_map_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_prok)
			print(run_cmd(cmd_run_dart_map_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))



			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_se_prok)
			print(run_cmd(cmd_run_dart_map_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print(run_cmd(cmd_run_dart_sort_bam))


			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print(run_cmd(cmd_run_dart_remove_unsorted_bam))

		# Run HISAT2
		####
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_euk)
			print(run_cmd(cmd_run_hisat2_pe_euk))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_euk)
			print(run_cmd(cmd_run_hisat2_se_euk))


		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_prok)
			print(run_cmd(cmd_run_hisat2_pe_euk))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_prok)
			print(run_cmd(cmd_run_hisat2_se_prok))


		# Run STAR

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_prok)
			print(run_cmd(cmd_run_star_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_prok)
			print(run_cmd(cmd_run_star_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))



		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_euk)
			print(run_cmd(cmd_run_star_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_euk)
			print(run_cmd(cmd_run_star_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print(run_cmd(cmd_star_bam_rename))

'''

	#####################################################################################################
class alingnReadSetsToGGAT(luigi.Task):

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
			return [alingnReadToGGAT(pre_process_reads=self.pre_process_reads,
								   annotation_file_type=self.annotation_file_type,
								   rnaseq_aligner=self.rnaseq_aligner,
								   sampleName=i)
						for i in [line.strip()
								  for line in
								  open(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))]]

		if self.read_library_type == "se":
			return [alingnReadToGGAT(pre_process_reads=self.pre_process_reads,
								   annotation_file_type=self.annotation_file_type,
								   rnaseq_aligner=self.rnaseq_aligner,
								   sampleName=i)
						for i in [line.strip()
								  for line in
								  open(os.path.join(os.getcwd(), "sample_list", "se_samples.lst"))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(), "task_logs", 'task.ref.based.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
				outfile.write('Reference based assembly finished at {t}'.format(t=timestamp))
'''

class clusterGGAT(luigi.Task):
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
	rnaseq_aligner = luigi.ChoiceParameter(choices=["star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF"], var_type=str)

	def requires(self):
		if self.read_library_type == "pe":
			return [alingnReadToGGAT(pre_process_reads=self.pre_process_reads,
									annotation_file_type=self.annotation_file_type,
									rnaseq_aligner=self.rnaseq_aligner,
									sampleName=i)
					for i in [line.strip()
							  for line in
							  open(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))]]

		if self.read_library_type == "se":
			return [alingnReadToGGAT(pre_process_reads=self.pre_process_reads,
									annotation_file_type=self.annotation_file_type,
									rnaseq_aligner=self.rnaseq_aligner,
									sampleName=i)
					for i in [line.strip()
							  for line in
							  open(os.path.join(os.getcwd(), "sample_list", "se_samples.lst"))]]

	def output(self):

		corset_folder = os.path.join(os.getcwd(), self.project_name,
									 "genome_guided_assembly_based_dea",
									 "ReadQuant",
									 "gg_trinity_" + self.rnaseq_aligner + "_corset_" +	 self.read_library_type)

		return {'out1': luigi.LocalTarget(corset_folder + "clusters.txt"),
				'out2': luigi.LocalTarget(corset_folder + "counts.txt")
				}

	def run(self):
		corset_folder = os.path.join(os.getcwd(), self.project_name,
									 "genome_guided_assembly_based_dea",
									 "ReadQuant",
									 "gg_trinity_" + self.rnaseq_aligner + "_corset_" + self.read_library_type)

		gg_transcript_map_folder = os.path.join(os.getcwd(), self.project_name,
								 "genome_guided_assembly_based_dea",
								 "rnaseq_transcript_alignment",
								 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

		gg_assembled_transcript_folder = os.path.join(os.getcwd(), self.project_name, "genome_guided_assembly",
													  "trinity_" + self.rnaseq_aligner + "_" + self.read_library_type + "/")

		super_transcript_folder = os.path.join(os.getcwd(), self.project_name, "genome_guided_assembly",
													  "trinity_" + self.rnaseq_aligner + "_" + self.read_library_type, "SuperTranscript" + "/")

		# Paired end reads

		input_group_file = os.path.join(os.getcwd(), "sample_list", "group.tsv")

		# Command to generate Rockhopper Paired-end input
		cmd_corset_read_input = prepare_corset_input(input_group_file)

		cmd_run_corset = "[ -d  {corset_folder} ] || mkdir -p {corset_folder}; cd {corset_folder}; " \
						 "corset -D 99999999999  {corset_read_input} {gg_transcript_map_folder}*.bam " \
			.format(corset_folder=corset_folder,
					corset_read_input=cmd_corset_read_input,
					gg_transcript_map_folder=gg_transcript_map_folder)

		cmd_run_supertrans = "python $(which Lace.py) " \
							 "--cores 1 " \
							 "{gg_assembled_transcript_folder}Trinity-GG.fasta " \
							 "{corset_folder}clusters.txt " \
							 "--outputDir {super_transcript_folder} " \
			.format(super_transcript_folder=super_transcript_folder, corset_folder=corset_folder,
					gg_assembled_transcript_folder=gg_assembled_transcript_folder)

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset)
		print(run_cmd(cmd_run_corset))

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_supertrans)
		print(run_cmd(cmd_run_supertrans))


class genomeGuidedTransAssemblyDEA(luigi.Task):
	adapter = GlobalParameter().adapter
	organism_domain = GlobalParameter().organism_domain
	read_library_type = GlobalParameter().read_library_type
	result_tag = GlobalParameter().result_tag

	project_name = luigi.Parameter(default="RNASeqAnalysis")
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF"], var_type=str)
	#rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	rnaseq_aligner = luigi.ChoiceParameter(choices=["star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	dea_method = luigi.ChoiceParameter(choices=["deseq2", "edger"], var_type=str)
	report_name = luigi.Parameter(default="Corset_DESeq2_HTML_Report")
	factor_of_intrest = luigi.Parameter(default="condititions",description="factor of intrest column of the target file (string [=condititions]). ")
	reference_condition = luigi.Parameter(default="control", description="reference biological condition.  (string [=control]")
	#target_file = luigi.Parameter(description="path to the design/target file. (string [=target.tsv]")
	alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
	p_adjust_method = luigi.ChoiceParameter(default="BH", description="p-value adjustment method.",choices=["BH", "BY"],var_type=str)
	fit_type = luigi.ChoiceParameter(default="parametric", description="mean-variance relationship.",  choices=["parametric", "local","mean"],var_type=str)
	size_factor = luigi.ChoiceParameter(default="median", description="method to estimate the size factors.", choices=["median", "shorth"],var_type=str)


	def requires(self):
		return [clusterGGAT(rnaseq_aligner=self.rnaseq_aligner,
							project_name=self.project_name,
							annotation_file_type=self.annotation_file_type,
						   pre_process_reads=self.pre_process_reads)]

	def output(self):

		#edgeResultFolter = GlobalParameter().basefolder + "/" + self.projectName + "_DEAnalysis/" + self.rnaseq_assembler + "_" +  self.sampleGroupFile + "_" + self.readType + "/" + "edgeR/"
		resultFolder = os.path.join(os.getcwd(), self.project_name,
									"genome_guided_assembly_based_dea",
									"DEAnalysis",
									self.dea_method + "_" + self.rnaseq_aligner + "_" +  self.result_tag + "_" + self.read_library_type + "/")



		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.dea_method == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.dea_method == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.dea_method == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.dea_method == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.dea_method == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.dea_method == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.dea_method == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.dea_method == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}




	def run(self):
		resultFolder = os.path.join(os.getcwd(), self.project_name,
									"genome_guided_assembly_based_dea",
									"DEAnalysis",
									self.dea_method + "_" + self.rnaseq_aligner + "_" + self.result_tag + "_" + self.read_library_type + "/")


		corset_folder = os.path.join(os.getcwd(), self.project_name,
									 "genome_guided_assembly_based_dea",
									 "ReadQuant",
									 "gg_trinity_" + self.rnaseq_aligner + "_corset_" + self.read_library_type)

		target_file = os.path.join(os.getcwd(),"sample_list","target.tsv")
		rmd_DESeq2File = os.path.expanduser(os.path.join(('~'), 'scriptome','tasks','utility',"PlotDESEQ2.Rmd"))
		rmd_edgeRFile = os.path.expanduser(os.path.join(('~'), 'scriptome', 'tasks', 'utility', "PlotEDGER.Rmd"))

		cmd_run_corset_deseq = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
						  "cd {resultFolder};" \
						  "corset_DESeq2.r " \
						  "-t {target_file} " \
						  "-q {corset_folder}counts.txt " \
						   "-v {factor_of_intrest} " \
						   "-c {reference_condition} " \
						   "-f {fit_type} " \
						   "-a {alpha} " \
						   "-p {p_adjust_method} " \
						   "-l {size_factor} " \
							   "-T {rmd_DESeq2File}" \
						.format(resultFolder=resultFolder,
								target_file=target_file,
								corset_folder=corset_folder,
								factor_of_intrest=self.factor_of_intrest,
								reference_condition=self.reference_condition,
								fit_type=self.fit_type,
								alpha=self.alpha,
								p_adjust_method=self.p_adjust_method,
								size_factor=self.size_factor,
								rmd_DESeq2File=rmd_DESeq2File)

		cmd_run_corset_edger = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
							   "cd {resultFolder};" \
							   "corset_edgeR.r " \
							   "-t {target_file} " \
							   "-q {corset_folder}counts.txt " \
							   "-v {factor_of_intrest} " \
							   "-c {reference_condition} " \
							   "-f {fit_type} " \
							   "-a {alpha} " \
							   "-p {p_adjust_method} " \
							   "-l {size_factor} " \
							   "-T {rmd_edgeRFile}" \
			.format(resultFolder=resultFolder,
					target_file=target_file,
					corset_folder=corset_folder,
					factor_of_intrest=self.factor_of_intrest,
					reference_condition=self.reference_condition,
					fit_type=self.fit_type,
					alpha=self.alpha,
					p_adjust_method=self.p_adjust_method,
					size_factor=self.size_factor,
					rmd_edgeRFile=rmd_edgeRFile)

		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.dea_method == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

##############################################################################################################
#Run EDGER
		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.dea_method == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

