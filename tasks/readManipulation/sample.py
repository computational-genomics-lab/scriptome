import luigi
import time
import os
import subprocess


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

def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)

createFolder("task_logs")

class byNumber(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["lr","pe", "mp","pe-mp","ilpe","ilmp","ilpe-ilmp"], var_type=str)

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
	number = luigi.IntParameter(description="number in terms of Kb of reads in sample data . (intiger)")


	def output(self):

		mp_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","mate-pair")
		pe_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","pair-end")
		pe_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","pair-end")
		mp_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","mate-pair")
		lr_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","long-read")

		if self.libType == "pe":
			return {'out1': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "mp":
			return {'out1': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "pe,mp":
			return {'out1': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilpe":
			return {'out1': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilmp":
			return {'out1': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilpe,ilmp":
			return {'out1': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R2.fastq")}

	def run(self):

		mp_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","mate-pair")
		pe_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","pair-end")
		pe_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","pair-end")
		mp_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","mate-pair")
		lr_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByNUM","long-read")

		sample_log_folder = os.path.join(os.getcwd(), self.projectName, "log","SampleReads_KB" + "/")


		cmd_smaple_pe = "[ -d  {pe_sample_read_folder} ] || mkdir -p {pe_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						"cd {pe_sample_read_folder}; reformat.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={paired_end_read_dir}{sampleName}_R1.{paired_end_read_suffix} " \
						"in2={paired_end_read_dir}{sampleName}_R2.{paired_end_read_suffix} " \
						"out={pe_sample_read_folder}{sampleName}_R1.fastq " \
						"out2={pe_sample_read_folder}{sampleName}_R2.fastq " \
						"samplereadstarget={number}k " \
						"ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						.format(sampleName=self.sampleName,
						paired_end_read_suffix=GlobalParameter().paired_end_read_suffix,
						pe_sample_read_folder=pe_sample_read_folder,
						sample_log_folder=sample_log_folder,
						Xmx=GlobalParameter().maxMemory,
						threads=GlobalParameter().threads,
						number=self.number,
						paired_end_read_dir=GlobalParameter().paired_end_read_dir)

		cmd_smaple_mp = "[ -d  {mp_sample_read_folder} ] || mkdir -p {mp_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						"cd {mp_sample_read_folder}; reformat.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={mate_pair_read_dir}{sampleName}_R1.{mate_pair_read_suffix} " \
						"in2={mate_pair_read_dir}{sampleName}_R2.{mate_pair_read_suffix} " \
						"out={mp_sample_read_folder}{sampleName}_R1.fastq " \
						"out2={mp_sample_read_folder}{sampleName}_R2.fastq " \
						"samplereadstarget={number}k " \
						"ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						.format(sampleName=self.sampleName,
						mate_pair_read_suffix=GlobalParameter().mate_pair_read_suffix,
						mp_sample_read_folder=mp_sample_read_folder,
						sample_log_folder=sample_log_folder,
						Xmx=GlobalParameter().maxMemory,
						threads=GlobalParameter().threads,
						number=self.number,
						mate_pair_read_dir=GlobalParameter().mate_pair_read_dir)

		cmd_sample_ilpe = "[ -d  {pe_il_sample_read_folder} ] || mkdir -p {pe_il_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						   "cd {pe_il_sample_read_folder}; reformat.sh " \
						  "-Xmx{Xmx}g " \
						  "threads={threads} " \
						  "in={paired_end_interleaved_read_dir}{sampleName}.{paired_end_interleaved_read_suffix} " \
						  "out={pe_il_sample_read_folder}{sampleName}_R1.fastq " \
						  "out2={pe_il_sample_read_folder}{sampleName}_R2.fastq " \
						  "samplereadstarget={number}k " \
						  "ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						  .format(sampleName=self.sampleName,
						  paired_end_interleaved_read_suffix=GlobalParameter().paired_end_interleaved_read_suffix,
						  pe_il_sample_read_folder=pe_il_sample_read_folder,
						  sample_log_folder=sample_log_folder,
						  Xmx=GlobalParameter().maxMemory,
						  threads=GlobalParameter().threads,
						  number=self.number,
						  paired_end_interleaved_read_dir=os.path.join(GlobalParameter().paired_end_interleaved_read_dir + "/"))

		cmd_sample_ilmp = "[ -d  {mp_il_sample_read_folder} ] || mkdir -p {mp_il_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						   "cd {mp_il_sample_read_folder}; reformat.sh " \
						   "-Xmx{Xmx}g " \
						   "threads={threads} " \
						   "in={mate_pair_interleaved_read_dir}{sampleName}.{mate_pair_interleaved_read_suffix} " \
						   "out={mp_il_sample_read_folder}{sampleName}_R1.fastq " \
						   "out2={mp_il_sample_read_folder}{sampleName}_R2.fastq " \
						   "ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						   "samplereadstarget={number}k " \
						   .format(sampleName=self.sampleName,
						   mate_pair_interleaved_read_suffix=GlobalParameter().mate_pair_interleaved_read_suffix,
						   mp_il_sample_read_folder=mp_il_sample_read_folder,
						   sample_log_folder=sample_log_folder,
						   Xmx=GlobalParameter().maxMemory,
						   threads=GlobalParameter().threads,
						   number=self.number,
						   mate_pair_interleaved_read_dir=os.path.join(GlobalParameter().mate_pair_interleaved_read_dir + "/"))


		if self.libType == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_pe)
			print(run_cmd(cmd_smaple_pe))

		if self.libType == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_mp)
			print(run_cmd(cmd_smaple_mp))

		if self.libType == "pe-mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_pe)
			print(run_cmd(cmd_smaple_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_mp)
			print(run_cmd(cmd_smaple_mp))

		if self.libType == "ilpe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilpe)
			print(run_cmd(cmd_sample_ilpe))

		if self.libType == "ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilmp)
			print(run_cmd(cmd_sample_ilmp))

		if self.libType == "ilpe-ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilpe)
			print(run_cmd(cmd_sample_ilpe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilmp)
			print(run_cmd(cmd_sample_ilmp))

class sampleByNumber(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	number = luigi.IntParameter(description="number in terms of Kb of reads in sample data . (intiger)")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["lr","pe", "mp", "pe-mp", "ilpe", "ilmp", "ilpe-ilmp"], var_type=str)
	def requires(self):
		if self.libType == "pe":
			return [byNumber(libType=self.libType,number=self.number,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), self.projectName, "pe_samples.lst")))]]
		if self.libType == "mp":
			return [byNumber(libType=self.libType,number=self.number,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), self.projectName, "mp_samples.lst")))]]
		if self.libType == "pe-mp":
			return [
						[byNumber(libType="pe",number=self.number,
							sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), self.projectName, "pe_samples.lst")))]],

						[byNumber(libType="mp",number=self.number,
							sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), self.projectName, "mp_samples.lst")))]]
				  ]


		if self.libType == "ilpe":
			return [byNumber(libType=self.libType,number=self.number,sampleName=i)
							for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), self.projectName, "pe_il_samples.lst")))]]
		if self.libType == "ilmp":
			return [byNumber(libType=self.libType,number=self.number,sampleName=i)
							for i in [line.strip()
							  for line in
								open((os.path.join(os.getcwd(), self.projectName, "mp_il_samples.lst")))]]
		if self.libType == "ilpe-ilmp":
			return [
						[byNumber(libType="ilpe",number=self.number,sampleName=i)
								for i in [line.strip()
									for line in
										open((os.path.join(os.getcwd(), self.projectName, "pe_il_samples.lst")))]],

						[byNumber(libType="ilmp",number=self.number, sampleName=i)
							for i in [line.strip()
								for line in
									open((os.path.join(os.getcwd(), self.projectName, "mp_il_samples.lst")))]]
				   ]


	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('data.downlaod.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('download finished at {t}'.format(t=timestamp))

#####################################################################################################
#BY FRACTION
#####################################################################################################

class byFraction(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["lr","pe", "mp","pe-mp","ilpe","ilmp","ilpe-ilmp"], var_type=str)

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
	fraction = luigi.FloatParameter(description="fraction of reads in sample data. Must be betwween 0.0 to 1.0 . (Float)")


	def output(self):

		mp_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","mate-pair")
		pe_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","pair-end")
		pe_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","pair-end")
		mp_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","mate-pair")
		lr_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","long-read")

		if self.libType == "pe":
			return {'out1': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "mp":
			return {'out1': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "pe,mp":
			return {'out1': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_sample_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilpe":
			return {'out1': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilmp":
			return {'out1': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "ilpe,ilmp":
			return {'out1': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_il_sample_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_il_sample_read_folder + self.sampleName + "_R2.fastq")}

	def run(self):

		mp_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","mate-pair" + "/")
		pe_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","pair-end" + "/")
		pe_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","pair-end" + "/")
		mp_il_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","mate-pair"+ "/")
		lr_sample_read_folder = os.path.join(os.getcwd(), self.projectName, "SampledReadsByFRAC","long-read"+ "/")

		sample_log_folder = os.path.join(os.getcwd(), self.projectName, "log","SampleReads_FRAC" + "/")


		cmd_smaple_pe = "[ -d  {pe_sample_read_folder} ] || mkdir -p {pe_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						"cd {pe_sample_read_folder}; reformat.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={paired_end_read_dir}{sampleName}_R1.{paired_end_read_suffix} " \
						"in2={paired_end_read_dir}{sampleName}_R2.{paired_end_read_suffix} " \
						"out={pe_sample_read_folder}{sampleName}_R1.fastq " \
						"out2={pe_sample_read_folder}{sampleName}_R2.fastq " \
						"samplerate={fraction} " \
						"ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						.format(sampleName=self.sampleName,
						paired_end_read_suffix=GlobalParameter().paired_end_read_suffix,
						pe_sample_read_folder=pe_sample_read_folder,
						sample_log_folder=sample_log_folder,
						Xmx=GlobalParameter().maxMemory,
						threads=GlobalParameter().threads,
						fraction=self.fraction,
						paired_end_read_dir=os.path.join(GlobalParameter().paired_end_read_dir + "/"))

		cmd_smaple_mp = "[ -d  {mp_sample_read_folder} ] || mkdir -p {mp_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						"cd {mp_sample_read_folder}; reformat.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={mate_pair_read_dir}{sampleName}_R1.{paired_end_interleaved_read_suffix} " \
						"in2={mate_pair_read_dir}{sampleName}_R2.{paired_end_interleaved_read_suffix} " \
						"out={mp_sample_read_folder}{sampleName}_R1.fastq " \
						"out2={mp_sample_read_folder}{sampleName}_R2.fastq " \
						"samplerate={fraction} " \
						"ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						.format(sampleName=self.sampleName,
						paired_end_interleaved_read_suffix=GlobalParameter().paired_end_interleaved_read_suffix,
						mp_sample_read_folder=mp_sample_read_folder,
						sample_log_folder=sample_log_folder,
						Xmx=GlobalParameter().maxMemory,
						threads=GlobalParameter().threads,
						fraction=self.fraction,
						mate_pair_read_dir=os.path.join(GlobalParameter().mate_pair_read_dir + "/"))


		cmd_sample_ilpe = "[ -d  {pe_il_sample_read_folder} ] || mkdir -p {pe_il_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						   "cd {pe_il_sample_read_folder}; reformat.sh " \
						  "-Xmx{Xmx}g " \
						  "threads={threads} " \
						  "in={paired_end_interleaved_read_dir}{sampleName}.{paired_end_read_suffix} " \
						  "out={pe_il_sample_read_folder}{sampleName}_R1.fastq " \
						  "out2={pe_il_sample_read_folder}{sampleName}_R2.fastq " \
						  "samplerate={fraction} " \
						  "ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						  .format(sampleName=self.sampleName,
						  paired_end_read_suffix=GlobalParameter().paired_end_read_suffix,
						  pe_il_sample_read_folder=pe_il_sample_read_folder,
						  sample_log_folder=sample_log_folder,
						  Xmx=GlobalParameter().maxMemory,
						  threads=GlobalParameter().threads,
						  fraction=self.fraction,
						  paired_end_interleaved_read_dir=os.path.join(GlobalParameter().paired_end_interleaved_read_dir + "/"))

		cmd_sample_ilmp = "[ -d  {mp_il_sample_read_folder} ] || mkdir -p {mp_il_sample_read_folder}; mkdir -p {sample_log_folder}; " \
						   "cd {mp_il_sample_read_folder}; reformat.sh " \
						   "-Xmx{Xmx}g " \
						   "threads={threads} " \
						   "in={mate_pair_interleaved_read_dir}{sampleName}.{mate_pair_interleaved_read_suffix} " \
						   "out={mp_il_sample_read_folder}{sampleName}_R1.fastq " \
						   "out2={mp_il_sample_read_folder}{sampleName}_R2.fastq " \
						   "ziplevel=9 2>{sample_log_folder}{sampleName}_sampling.log" \
						   "samplerate={fraction} " \
						   .format(sampleName=self.sampleName,
						   mate_pair_interleaved_read_suffix=GlobalParameter().mate_pair_interleaved_read_suffix,
						   mp_il_sample_read_folder=mp_il_sample_read_folder,
						   sample_log_folder=sample_log_folder,
						   Xmx=GlobalParameter().maxMemory,
						   threads=GlobalParameter().threads,
						   fraction=self.fraction,
						   mate_pair_interleaved_read_dir=os.path.join(GlobalParameter().mate_pair_interleaved_read_dir + "/"))


		if self.libType == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_pe)
			print(run_cmd(cmd_smaple_pe))

		if self.libType == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_mp)
			print(run_cmd(cmd_smaple_mp))

		if self.libType == "pe-mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_pe)
			print(run_cmd(cmd_smaple_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_smaple_mp)
			print(run_cmd(cmd_smaple_mp))

		if self.libType == "ilpe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilpe)
			print(run_cmd(cmd_sample_ilpe))

		if self.libType == "ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilmp)
			print(run_cmd(cmd_sample_ilmp))

		if self.libType == "ilpe-ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilpe)
			print(run_cmd(cmd_sample_ilpe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_sample_ilmp)
			print(run_cmd(cmd_sample_ilmp))

class sampleByFraction(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	fraction = luigi.FloatParameter(description="fraction of reads in sample data. Must be betwween 0.0 to 1.0 . (Float)")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["lr","pe", "mp", "pe-mp", "ilpe", "ilmp", "ilpe-ilmp"], var_type=str)
	def requires(self):
		if self.libType == "pe":
			return [byFraction(libType=self.libType,fraction=self.fraction,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		if self.libType == "mp":
			return [byFraction(libType=self.libType,fraction=self.fraction,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
		if self.libType == "pe-mp":
			return [
						[byFraction(libType="pe",fraction=self.fraction,
							sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

						[byFraction(libType="mp",fraction=self.fraction,
							sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
				  ]


		if self.libType == "ilpe":
			return [byFraction(libType=self.libType,fraction=self.fraction,sampleName=i)
							for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]]
		if self.libType == "ilmp":
			return [byFraction(libType=self.libType,fraction=self.fraction,sampleName=i)
							for i in [line.strip()
							  for line in
								open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]
		if self.libType == "ilpe-ilmp":
			return [
						[byFraction(libType="ilpe",fraction=self.fraction,sampleName=i)
								for i in [line.strip()
									for line in
										open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]],

						[byFraction(libType="ilmp",fraction=self.fraction,sampleName=i)
							for i in [line.strip()
								for line in
									open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]
				   ]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(), "task_logs", 'task.data.sampling.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('data sampling finished at {t}'.format(t=timestamp))