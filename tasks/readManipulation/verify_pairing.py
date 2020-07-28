import luigi
import time
import os
import subprocess


class GlobalParameter(luigi.Config):
	pairedendDir = luigi.Parameter()
	matepairDir = luigi.Parameter()
	peInterleavedDir = luigi.Parameter()
	mpInterleavedDir = luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	suffix = luigi.Parameter()

def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output


class verifyPair(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["pe", "mp","pe-mp","ilpe","ilmp","ilpe-ilmp"], var_type=str)

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")

	def output(self):

		verify_pairing_log_folder = os.path.join(os.getcwd(), "log","verify_pairing" + "/")



		if self.libType == "pe":
			return {'out1': luigi.LocalTarget(verify_pairing_log_folder + self.sampleName + "_verify_pairing.log")}

		if self.libType == "mp":
			return {'out1': luigi.LocalTarget(verify_pairing_log_folder + self.sampleName + "_verify_pairing.log")}

		if self.libType == "pe,mp":
			return {'out1': luigi.LocalTarget(verify_pairing_log_folder + self.sampleName + "_verify_pairing.log")}


		if self.libType == "ilpe":
			return {'out1': luigi.LocalTarget(verify_pairing_log_folder + self.sampleName + "_verify_pairing.log")}

		if self.libType == "ilmp":
			return {'out1': luigi.LocalTarget(verify_pairing_log_folder + self.sampleName + "_verify_pairing.log")}

		if self.libType == "ilpe,ilmp":
			return {'out1': luigi.LocalTarget(verify_pairing_log_folder + self.sampleName + "_verify_pairing.log")}

	def run(self):

		verify_pairing_log_folder = os.path.join(os.getcwd(), "log","verify_pairing" + "/")

		cmd_verify_pe = "[ -d  {verify_pairing_log_folder} ] || mkdir -p {verify_pairing_log_folder}; " \
						"cd {verify_pairing_log_folder}; reformat.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={pairedendDir}{sampleName}_R1.{suffix} " \
						"in2={pairedendDir}{sampleName}_R2.{suffix} " \
						"vpair 2>{verify_pairing_log_folder}{sampleName}_verify_pairing.log " \
			.format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					verify_pairing_log_folder=verify_pairing_log_folder,
					pairedendDir=os.path.join(GlobalParameter().pairedendDir + "/"))

		cmd_verify_mp = "[ -d  {verify_pairing_log_folder} ] || mkdir -p {verify_pairing_log_folder}; " \
						"cd {verify_pairing_log_folder}; reformat.sh " \
						"-Xmx{Xmx}g " \
						"threads={threads} " \
						"in1={matepairDir}{sampleName}_R1.{suffix} " \
						"in2={matepairDir}{sampleName}_R2.{suffix} " \
						"vpair 2>{verify_pairing_log_folder}{sampleName}_verify_pairing.log " \
						.format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					verify_pairing_log_folder=verify_pairing_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					matepairDir=os.path.join(GlobalParameter().matepairDir + "/"))

		cmd_verify_pe_il = "[ -d  {verify_pairing_log_folder} ] || mkdir -p {verify_pairing_log_folder}; " \
						   "cd {verify_pairing_log_folder}; reformat.sh " \
						  "-Xmx{Xmx}g " \
						  "threads={threads} " \
						  "in={peInterleavedDir}{sampleName}.{suffix} " \
						  "verifypairing 2>{verify_pairing_log_folder}{sampleName}_verify_pairing.log " \
						  .format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					verify_pairing_log_folder=verify_pairing_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					peInterleavedDir=os.path.join(GlobalParameter().peInterleavedDir + "/"))

		cmd_verify_mp_il = "[ -d  {verify_pairing_log_folder} ] || mkdir -p {verify_pairing_log_folder}; " \
						   "cd {verify_pairing_log_folder}; reformat.sh " \
						   "-Xmx{Xmx}g " \
						   "threads={threads} " \
						   "in={mpInterleavedDir}{sampleName}.{suffix} " \
						   "verifypairing 2>{verify_pairing_log_folder}{sampleName}_verify_pairing.log " \
						   .format(sampleName=self.sampleName,
					suffix=GlobalParameter().suffix,
					verify_pairing_log_folder=verify_pairing_log_folder,
					Xmx=GlobalParameter().maxMemory,
					threads=GlobalParameter().threads,
					mpInterleavedDir=os.path.join(GlobalParameter().mpInterleavedDir + "/"))


		if self.libType == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe)
			print(run_cmd(cmd_verify_pe))

		if self.libType == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_mp)
			print(run_cmd(cmd_verify_mp))

		if self.libType == "pe-mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe)
			print(run_cmd(cmd_verify_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_mp)
			print(run_cmd(cmd_verify_mp))

		if self.libType == "ilpe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe_il)
			print(run_cmd(cmd_verify_pe_il))

		if self.libType == "ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_mp_il)
			print(run_cmd(cmd_verify_mp_il))

		if self.libType == "ilpe-ilmp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe_il)
			print(run_cmd(cmd_verify_pe_il))
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_mp_il)
			print(run_cmd(cmd_verify_mp_il))



class verifyPairSamples(luigi.Task):
	projectName = luigi.Parameter(default="ReadManipulation")
	libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
												"'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
												"'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
									choices=["pe", "mp", "pe-mp", "ilpe", "ilmp", "ilpe-ilmp"], var_type=str)


	def requires(self):
		if self.libType == "pe":
			return [verifyPair(libType=self.libType,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		if self.libType == "mp":
			return [verifyPair(libType=self.libType,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(),"sample_list", "mp_samples.lst")))]]
		if self.libType == "pe-mp":
			return [
						[verifyPair(libType="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[verifyPair(libType="mp", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","mp_samples.lst")))]]
				  ]

		if self.libType == "ilpe":
			return [verifyPair(libType=self.libType,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]]
		if self.libType == "ilmp":
			return [verifyPair(libType=self.libType,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]
		if self.libType == "ilpe-ilmp":
			return [
				[verifyPair(libType="ilpe", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(),"sample_list",  "pe_samples.lst")))]],

				[verifyPair(libType="ilmp", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(),"sample_list",  "mp_samples.lst")))]]
			]


	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('task.verifypairing.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('pair verification finished at {t}'.format(t=timestamp))
