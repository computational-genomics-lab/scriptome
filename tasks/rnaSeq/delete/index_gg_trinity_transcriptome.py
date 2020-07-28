import luigi
import os
import subprocess
from Bio import SeqIO
import math
from tasks.assembly.gg_trinity import ggTrinity

def run_cmd(cmd):
    p = subprocess.Popen(cmd, bufsize=-1,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
    output = p.communicate()[0]
    return output


class GlobalParameter(luigi.Config):
	genome_suffix=luigi.Parameter()
	read_library_type=luigi.Parameter()
	organism_domain=luigi.Parameter()
	genome_name=luigi.Parameter()
	genome_dir=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()

def genomeSAindexNbases(genome_fasta):
    with open(genome_fasta) as genome:
        totalsize=0
        for rec in SeqIO.parse(genome, 'fasta'):
            totalsize = totalsize + len(rec)
        log2 = math.log(totalsize, 2.0)
        index = (log2/2) - 1
        gsanb = int(index)
        return gsanb

class indexTrinityTranscript(luigi.Task):
	project_name = luigi.Parameter(default="RNASeqAnalysis")
	genome_name = GlobalParameter().genome_name
	genome_suffix= GlobalParameter().genome_suffix
	genome_dir=GlobalParameter().genome_dir
	organism_domain = GlobalParameter().organism_domain
	read_library_type= GlobalParameter().read_library_type
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star","hisat2","dart", "segemehl","bowtie2"],var_type=str)
	annotation_file_type = luigi.ChoiceParameter(choices=["GFF","GTF"],var_type=str)
	pre_process_reads=luigi.ChoiceParameter(choices=["yes","no"],var_type=str)

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
		transcriptIndexFolder=os.path.join(os.getcwd(), self.project_name, "transcript_index", self.genome_name +"_"+ self.rnaseq_aligner +"_index" +
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

		transcriptomeFolder=os.path.join(os.getcwd(), self.project_name, "genomeguided_assembly",
													  "trinity_" + self.rnaseq_aligner + "_" + self.read_library_type + "/")

		transcriptIndexFolder = os.path.join(os.getcwd(), self.project_name, "transcript_index",
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