import luigi
import time
import os
import subprocess

def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)

createFolder("task_logs")

class GlobalParameter(luigi.Config):
	paired_end_read_dir = luigi.Parameter()
	mate_pair_read_dir = luigi.Parameter()
	single_end_read_dir = luigi.Parameter()
	paired_end_interleaved_read_dir = luigi.Parameter()
	mate_pair_interleaved_read_dir = luigi.Parameter()
	long_read_dir = luigi.Parameter()
	trusted_contigs_dir = luigi.Parameter()
	genome_dir = luigi.Parameter()
	sequencing_read_type = luigi.Parameter()
	read_library_type = luigi.Parameter()
	genome_size = luigi.Parameter()
	organism_domain = luigi.Parameter()

	paired_end_read_suffix = luigi.Parameter()
	single_end_read_suffix = luigi.Parameter()
	paired_end_interleaved_read_suffix = luigi.Parameter()
	mate_pair_read_suffix = luigi.Parameter()
	mate_pair_interleaved_read_suffix = luigi.Parameter()
	long_read_suffix = luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()

def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output


class rePair(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["pe", "mp","pe-mp","ilpe","ilmp","ilpe-ilmp"], var_type=str)

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")

	def output(self):

		mp_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","mate-pair" + "/")
		pe_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","pair-end" + "/")
		pe_il_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","pair-end" + "/")
		mp_il_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","mate-pair" + "/")

		if self.libType == "pe":
			return {'out1': luigi.LocalTarget(pe_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_repair_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "mp":
			return {'out1': luigi.LocalTarget(mp_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_repair_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "pe,mp":
			return {'out1': luigi.LocalTarget(pe_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_repair_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_repair_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilpe":
			return {'out1': luigi.LocalTarget(pe_il_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_il_repair_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilmp":
			return {'out1': luigi.LocalTarget(mp_il_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_il_repair_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilpe,ilmp":
			return {'out1': luigi.LocalTarget(pe_il_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_il_repair_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_il_repair_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_il_repair_read_folder + self.sampleName + "_R2.fastq")}

	def run(self):

		mp_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","mate-pair" + "/")
		pe_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","pair-end" + "/")
		pe_il_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","pair-end" + "/")
		mp_il_repair_read_folder = os.path.join(os.getcwd(), self.projectName, "RePairedReads","mate-pair" + "/")
		re_pair_log_folder = os.path.join(os.getcwd(), self.projectName, "log","re_pair_reads" + "/")


		cmd_repair_pe = "[ -d  {pe_repair_read_folder} ] || mkdir -p {pe_repair_read_folder}; mkdir -p {re_pair_log_folder}; " \
						"cd {pe_repair_read_folder}; repair.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={pairedendDir}{sampleName}_R1.{suffix} " \
						"in2={pairedendDir}{sampleName}_R2.{suffix} " \
						"out={pe_repair_read_folder}{sampleName}_R1.fastq " \
						"out2={pe_repair_read_folder}{sampleName}_R2.fastq " \
						"outsingle={pe_repair_read_folder}{sampleName}_singleton.fastq " \
						"ziplevel=9 2>{re_pair_log_folder}{sampleName}_re_pair.log" \
			.format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					pe_repair_read_folder=pe_repair_read_folder,
					re_pair_log_folder=re_pair_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					pairedendDir=os.path.join(GlobalParameter().pairedendDir + "/"))

		cmd_repair_mp = "[ -d  {mp_repair_read_folder} ] || mkdir -p {mp_repair_read_folder}; mkdir -p {re_pair_log_folder}; " \
						"cd {mp_repair_read_folder}; repair.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={matepairDir}{sampleName}_R1.{suffix} " \
						"in2={matepairDir}{sampleName}_R2.{suffix} " \
						"out={mp_repair_read_folder}{sampleName}_R1.fastq " \
						"out2={mp_repair_read_folder}{sampleName}_R2.fastq " \
						"outsingle={mp_repair_read_folder}{sampleName}_singleton.fastq " \
						"ziplevel=9 2>{re_pair_log_folder}{sampleName}_re_pair.log" \
			.format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					mp_repair_read_folder=mp_repair_read_folder,
					re_pair_log_folder=re_pair_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					matepairDir=os.path.join(GlobalParameter().matepairDir + "/"))


		cmd_repair_pe_il = "[ -d  {pe_il_repair_read_folder} ] || mkdir -p {pe_il_repair_read_folder}; mkdir -p {re_pair_log_folder}; " \
						   "cd {pe_il_repair_read_folder}; repair.sh " \
						  "-Xmx{Xmx}g " \
						  "threads={threads} " \
						  "in={peInterleavedDir}{sampleName}.{suffix} " \
						  "out={pe_il_repair_read_folder}{sampleName}_R1.fastq " \
						  "out2={pe_il_repair_read_folder}{sampleName}_R2.fastq " \
						  "outsingle={pe_il_repair_read_folder}{sampleName}_singleton.fastq " \
						  "ziplevel=9 2>{re_pair_log_folder}{sampleName}_re_pair.log" \
			.format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					pe_il_repair_read_folder=pe_il_repair_read_folder,re_pair_log_folder=re_pair_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					peInterleavedDir=os.path.join(GlobalParameter().peInterleavedDir + "/"))

		cmd_repair_mp_il = "[ -d  {mp_il_repair_read_folder} ] || mkdir -p {mp_il_repair_read_folder}; mkdir -p {re_pair_log_folder}; " \
						   "cd {mp_il_repair_read_folder}; repair.sh " \
						   "-Xmx{Xmx}g " \
						   "threads={threads} " \
						   "in={mpInterleavedDir}{sampleName}.{suffix} " \
						   "out={mp_il_repair_read_folder}{sampleName}_R1.fastq " \
						   "out2={mp_il_repair_read_folder}{sampleName}_R2.fastq " \
						   "outsingle={mp_il_repair_read_folder}{sampleName}_singleton.fastq " \
						   "ziplevel=9 2>{re_pair_log_folder}{sampleName}_re_pair.log" \
			.format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					mp_il_repair_read_folder=mp_il_repair_read_folder,
					re_pair_log_folder=re_pair_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					mpInterleavedDir=os.path.join(GlobalParameter().mpInterleavedDir + "/"))


		if self.libType == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_pe)
			print(run_cmd(cmd_repair_pe))

		if self.libType == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_mp)
			print(run_cmd(cmd_repair_mp))

		if self.libType == "pe-mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_pe)
			print(run_cmd(cmd_repair_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_mp)
			print(run_cmd(cmd_repair_mp))

		if self.libType == "ilpe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_pe_il)
			print(run_cmd(cmd_repair_pe_il))

		if self.libType == "ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_mp_il)
			print(run_cmd(cmd_repair_mp_il))

		if self.libType == "ilpe-ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_pe_il)
			print(run_cmd(cmd_repair_pe_il))
			print("****** NOW RUNNING COMMAND ******: " + cmd_repair_mp_il)
			print(run_cmd(cmd_repair_mp_il))

class rePairSamples(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["pe", "mp", "pe-mp", "ilpe", "ilmp", "ilpe-ilmp"], var_type=str)


	def requires(self):
		if self.libType == "pe":
			return [rePair(libType=self.libType,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		if self.libType == "mp":
			return [rePair(libType=self.libType,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(),"sample_list", "mp_samples.lst")))]]
		if self.libType == "pe-mp":
			return [
						[rePair(libType="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

						[rePair(libType="mp", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
				  ]


		if self.libType == "ilpe":
			return [rePair(libType=self.libType,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]]
		if self.libType == "ilmp":
			return [rePair(libType=self.libType,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]
		if self.libType == "ilpe-ilmp":
			return [
				[rePair(libType="ilpe", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]],

				[rePair(libType="ilmp", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]
			]


	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('data.downlaod.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('download finished at {t}'.format(t=timestamp))