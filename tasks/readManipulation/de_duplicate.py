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


class deDup(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["Se","Pe", "Mp", "PeMp","PeMpSe", "Ilpe", "Ilmp", "IlpeIlmp", "IlpeIlmpSe"], var_type=str)

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")

	def output(self):

		mp_dedup_read_folder = os.path.join(os.getcwd(), self.projectName, "deDuplicatedReads","mate-pair" + "/")
		pe_dedup_read_folder = os.path.join(os.getcwd(), self.projectName, "deDuplicatedReads","pair-end" + "/")
		se_dedup_read_folder = os.path.join(os.getcwd(), self.projectName, "deDuplicatedReads","single-end" + "/")

		if self.libType == "Pe":
			return {'out1': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "Mp":
			return {'out1': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "Se":
			return {'out': luigi.LocalTarget(se_dedup_read_folder + self.sampleName + ".fastq")}

		if self.libType == "PeMp":
			return {'out1': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "PeMpSe":
			return {'out1': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R2.fastq"),
					'out5': luigi.LocalTarget(se_dedup_read_folder + self.sampleName + ".fastq")}

		if self.libType == "Ilpe":
			return {'out1': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "Ilmp":
			return {'out1': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R2.fastq")}

		if self.libType == "IlpeIlmp":
			return {'out1': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R2.fastq")}
		if self.libType == "IlpeIlmpSe":
			return {'out1': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_dedup_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_dedup_read_folder + self.sampleName + "_R2.fastq"),
					'out5': luigi.LocalTarget(se_dedup_read_folder + self.sampleName + ".fastq")}

	def run(self):

		mp_dedup_read_folder = os.path.join(os.getcwd(), self.projectName, "deDuplicatedReads","mate-pair" + "/")
		se_dedup_read_folder = os.path.join(os.getcwd(), self.projectName, "deDuplicatedReads","single-end" + "/")
		pe_dedup_read_folder = os.path.join(os.getcwd(), self.projectName, "deDuplicatedReads","pair-end" + "/")
		dedupicate_log_folder = os.path.join(os.getcwd(), self.projectName, "log","deduplicate_reads" + "/")

		cmd_dedup_pe = "[ -d  {pe_dedup_read_folder} ] || mkdir -p {pe_dedup_read_folder}; mkdir -p {dedupicate_log_folder}; " \
						"cd {pe_dedup_read_folder}; clumpify.sh dedupe=t " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={paired_end_read_dir}{sampleName}_R1.{paired_end_read_suffix} " \
						"in2={paired_end_read_dir}{sampleName}_R2.{paired_end_read_suffix} " \
						"out={pe_dedup_read_folder}{sampleName}_R1.fastq " \
						"out2={pe_dedup_read_folder}{sampleName}_R2.fastq " \
						"ziplevel=9 2>{dedupicate_log_folder}{sampleName}_dedup.log " \
						.format(sampleName=self.sampleName,
						paired_end_read_suffix=GlobalParameter().paired_end_read_suffix,
						pe_dedup_read_folder=pe_dedup_read_folder,dedupicate_log_folder=dedupicate_log_folder,
						Xmx=GlobalParameter().maxMemory,
						threads=GlobalParameter().threads,
						paired_end_read_dir=os.path.join(GlobalParameter().paired_end_read_dir + "/"))

		cmd_dedup_mp = "[ -d  {mp_dedup_read_folder} ] || mkdir -p {mp_dedup_read_folder}; mkdir -p {dedupicate_log_folder}; " \
						"cd {mp_dedup_read_folder}; clumpify.sh dedupe=t " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={mate_pair_read_dir}{sampleName}_R1.{mate_pair_read_suffix} " \
						"in2={mate_pair_read_dir}{sampleName}_R2.{mate_pair_read_suffix} " \
						"out={mp_dedup_read_folder}{sampleName}_R1.fastq " \
						"out2={mp_dedup_read_folder}{sampleName}_R2.fastq " \
						"ziplevel=9 2>{dedupicate_log_folder}{sampleName}_dedup.log " \
						.format(sampleName=self.sampleName,
					mate_pair_read_suffix=GlobalParameter().mate_pair_read_suffix,
					mp_dedup_read_folder=mp_dedup_read_folder,
					dedupicate_log_folder=dedupicate_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					mate_pair_read_dir=os.path.join(GlobalParameter().mate_pair_read_dir + "/"))

		cmd_dedup_se = "[ -d  {se_dedup_read_folder} ] || mkdir -p {se_dedup_read_folder}; mkdir -p {se_dedup_read_folder}; " \
						   "cd {se_dedup_read_folder}; clumpify.sh dedupe=t " \
						   "-Xmx{Xmx}g " \
						   "threads={threads} " \
						   "in={single_end_read_dir}{sampleName}.{single_end_read_suffix} " \
						   "out={se_dedup_read_folder}{sampleName}.fastq " \
						   "ziplevel=9 2>{dedupicate_log_folder}{sampleName}_dedup.log " \
						   .format(sampleName=self.sampleName,
					single_end_read_suffix=GlobalParameter().single_end_read_suffix,
					se_dedup_read_folder=se_dedup_read_folder,
					dedupicate_log_folder=dedupicate_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					single_end_read_dir=os.path.join(GlobalParameter().single_end_read_dir + "/"))

		cmd_dedup_pe_il = "[ -d  {pe_dedup_read_folder} ] || mkdir -p {pe_dedup_read_folder}; mkdir -p {dedupicate_log_folder}; " \
						   "cd {pe_dedup_read_folder}; clumpify.sh dedupe=t " \
						  "-Xmx{Xmx}g " \
						  "threads={threads} " \
						  "in={paired_end_interleaved_read_dir}{sampleName}.{paired_end_interleaved_read_suffix} " \
						  "out={pe_dedup_read_folder}{sampleName}_R1.fastq " \
						  "out2={pe_dedup_read_folder}{sampleName}_R2.fastq " \
						  "ziplevel=9 2>{dedupicate_log_folder}{sampleName}_dedup.log " \
						  .format(sampleName=self.sampleName,
					paired_end_interleaved_read_suffix=GlobalParameter().paired_end_interleaved_read_suffix,
					pe_dedup_read_folder=pe_dedup_read_folder,
					dedupicate_log_folder=dedupicate_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					paired_end_interleaved_read_dir=os.path.join(GlobalParameter().paired_end_interleaved_read_dir + "/"))

		cmd_dedup_mp_il = "[ -d  {mp_dedup_read_folder} ] || mkdir -p {mp_dedup_read_folder}; mkdir -p {dedupicate_log_folder}; " \
						   "cd {mp_dedup_read_folder}; clumpify.sh dedupe=t " \
						   "-Xmx{Xmx}g " \
						   "threads={threads} " \
						   "in={mate_pair_interleaved_read_dir}{sampleName}.{mate_pair_interleaved_read_suffix} " \
						   "out={mp_dedup_read_folder}{sampleName}_R1.fastq " \
						   "out2={mp_dedup_read_folder}{sampleName}_R2.fastq " \
						   "ziplevel=9 2>{dedupicate_log_folder}{sampleName}_dedup.log " \
						   .format(sampleName=self.sampleName,
					mate_pair_interleaved_read_suffix=GlobalParameter().mate_pair_interleaved_read_suffix,
					mp_dedup_read_folder=mp_dedup_read_folder,
					dedupicate_log_folder=dedupicate_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					mate_pair_interleaved_read_dir=os.path.join(GlobalParameter().mate_pair_interleaved_read_dir + "/"))


		if self.libType == "Pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_pe)
			print(run_cmd(cmd_dedup_pe))

		if self.libType == "Se":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_se)
			print(run_cmd(cmd_dedup_se))

		if self.libType == "Mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_mp)
			print(run_cmd(cmd_dedup_mp))

		if self.libType == "PeMp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_pe)
			print(run_cmd(cmd_dedup_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_mp)
			print(run_cmd(cmd_dedup_mp))

		if self.libType == "PeMpSe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_pe)
			print(run_cmd(cmd_dedup_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_mp)
			print(run_cmd(cmd_dedup_mp))
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_se)
			print(run_cmd(cmd_dedup_se))

		if self.libType == "Ilpe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_pe_il)
			print(run_cmd(cmd_dedup_pe_il))

		if self.libType == "Ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_mp_il)
			print(run_cmd(cmd_dedup_mp_il))

		if self.libType == "IlpeIlmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_pe_il)
			print(run_cmd(cmd_dedup_pe_il))
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_mp_il)
			print(run_cmd(cmd_dedup_mp_il))
		if self.libType == "IlpeIlmpSe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_pe_il)
			print(run_cmd(cmd_dedup_pe_il))
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_mp_il)
			print(run_cmd(cmd_dedup_mp_il))
			print("****** NOW RUNNING COMMAND ******: " + cmd_dedup_se)
			print(run_cmd(cmd_dedup_se))


class deDupSamples(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["Se","Pe", "Mp", "PeMp","PeMpSe", "Ilpe", "Ilmp", "IlpeIlmp", "IlpeIlmpSe"], var_type=str)


	def requires(self):

		if self.libType == "Pe":
			return [deDup(libType=self.libType,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		if self.libType == "Mp":
			return [deDup(libType=self.libType,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]

		if self.libType == "Se":
			return [deDup(libType=self.libType,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]


		if self.libType == "PeMp":
			return [
						[deDup(libType="Pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

						[deDup(libType="Mp", sampleName=i)
								for i in [line.strip()
										  for line in
											open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
				  ]

		if self.libType == "PeMpSe":
			return [
						[deDup(libType="Pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

						[deDup(libType="Mp", sampleName=i)
								for i in [line.strip()
										  for line in
											open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]],
						[deDup(libType="Se",
						   sampleName=i)
								for i in [line.strip()
											for line in
											open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]
				  ]

		if self.libType == "Ilpe":
			return [deDup(libType=self.libType,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]]


		if self.libType == "Ilmp":
			return [deDup(libType=self.libType,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]


		if self.libType == "IlpeIlmpSe":
			return [
					[deDup(libType="Ilpe", sampleName=i)
							for i in [line.strip()
								for line in
									open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]],

					[deDup(libType="Ilmp", sampleName=i)
							for i in [line.strip()
								for line in
									open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]],
					[deDup(libType="Se", sampleName=i)
							for i in [line.strip()
								for line in
									open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]
					]
	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.deduplication.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('read deduplication finished at {t}'.format(t=timestamp))