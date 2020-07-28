import luigi
import time
import os
import subprocess


class GlobalParameter(luigi.Config):
	projectName = luigi.Parameter(default="ReadManipulation")
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


class split(luigi.Task):
    projectName = luigi.Parameter(default="ReadManipulation")

    libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
                                                "'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
                                                "'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
                                    choices=["pe", "mp","pe-mp","ilpe","ilmp","ilpe-ilmp"], var_type=str)

    sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
    parts = luigi.IntParameter(description="Number of parts")

    def output(self):

        pe_split_read_folder = os.path.join(os.getcwd(), self.projectName, "SplittedReads","pair-end" + "/")
        mp_split_read_folder = os.path.join(os.getcwd(), self.projectName, "SplittedReads","mate-pair" + "/")


        if self.libType == "pe":
            return {'out1': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out2': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R2.fastq")}

        if self.libType == "mp":
            return {'out1': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out2': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R2.fastq")}

        if self.libType == "pe,mp":
            return {'out1': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out2': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R2.fastq"),
                    'out3': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out4': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R2.fastq")}

        if self.libType == "ilpe":
            return {'out1': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out2': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R2.fastq")}

        if self.libType == "ilmp":
            return {'out1': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out2': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R2.fastq")}

        if self.libType == "ilpe,ilmp":
            return {'out1': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out2': luigi.LocalTarget(pe_split_read_folder + self.sampleName + "_P1" + "_R2.fastq"),
                    'out3': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R1.fastq"),
                    'out4': luigi.LocalTarget(mp_split_read_folder + self.sampleName + "_P1" + "_R2.fastq")}
    def run(self):

        pe_split_read_folder = os.path.join(os.getcwd(), self.projectName, "SplittedReads","pair-end" + "/")
        mp_split_read_folder = os.path.join(os.getcwd(), self.projectName, "SplittedReads","mate-pair" + "/")
        split_log_folder = os.path.join(os.getcwd(), self.projectName, "log","split_reads" + "/")

        cmd_split_pe = "[ -d  {pe_split_read_folder} ] || mkdir -p {pe_split_read_folder}; mkdir -p {split_log_folder}; " \
                        "cd {pe_split_read_folder}; partition.sh " \
                        "-Xmx{Xmx}g " \
                        "threads={threads} " \
                        "in1={pairedendDir}{sampleName}_R1.{suffix} " \
                        "in2={pairedendDir}{sampleName}_R2.{suffix} " \
                        "out={pe_split_read_folder}{sampleName}_P%_R1.fastq " \
                        "out2={pe_split_read_folder}{sampleName}_P%_R2.fastq " \
                        "ways={parts} 2>{split_log_folder}{sampleName}_split.log " \
                        .format(sampleName=self.sampleName,
                        suffix=GlobalParameter().suffix,
                        pe_split_read_folder=pe_split_read_folder,split_log_folder=split_log_folder,
                        Xmx=GlobalParameter().maxMemory,
                        threads=GlobalParameter().threads,
                        parts=self.parts,
                        pairedendDir=os.path.join(GlobalParameter().pairedendDir + "/"))

        cmd_split_mp = "[ -d  {mp_split_read_folder} ] || mkdir -p {mp_split_read_folder}; mkdir -p {split_log_folder}; " \
                        "cd {mp_split_read_folder}; partition.sh " \
                        "-Xmx{Xmx}g " \
                        "threads={threads} " \
                        "in1={matepairDir}{sampleName}_R1.{suffix} " \
                        "in2={matepairDir}{sampleName}_R2.{suffix} " \
                        "out={mp_split_read_folder}{sampleName}_P%_R1.fastq " \
                        "out2={mp_split_read_folder}{sampleName}_P%_R2.fastq " \
                        "ways={parts} 2>{split_log_folder}{sampleName}_split.log " \
                        .format(sampleName=self.sampleName,
                    suffix=GlobalParameter().suffix,
                    mp_split_read_folder=mp_split_read_folder,
                    split_log_folder=split_log_folder,
                    Xmx=GlobalParameter().maxMemory,
                    threads=GlobalParameter().threads,
                    parts=self.parts,
                    matepairDir=os.path.join(GlobalParameter().matepairDir + "/"))

        cmd_split_pe_il = "[ -d  {pe_split_read_folder} ] || mkdir -p {pe_split_read_folder}; mkdir -p {split_log_folder}; " \
                           "cd {pe_split_read_folder}; partition.sh " \
                          "-Xmx{Xmx}g " \
                          "threads={threads} " \
                          "in={peInterleavedDir}{sampleName}.{suffix} " \
                          "out={pe_split_read_folder}{sampleName}_P%_R1.fastq " \
                          "out2={pe_split_read_folder}{sampleName}_P%_R2.fastq " \
                          "ways={parts} 2>{split_log_folder}{sampleName}_split.log " \
                          .format(sampleName=self.sampleName,
                    suffix=GlobalParameter().suffix,
                    pe_split_read_folder=pe_split_read_folder,
                    split_log_folder=split_log_folder,
                    Xmx=GlobalParameter().maxMemory,
                    threads=GlobalParameter().threads,
                    parts=self.parts,
                    peInterleavedDir=os.path.join(GlobalParameter().peInterleavedDir + "/"))

        cmd_split_mp_il = "[ -d  {mp_split_read_folder} ] || mkdir -p {mp_split_read_folder}; mkdir -p {split_log_folder}; " \
                           "cd {mp_split_read_folder}; partition.sh " \
                           "-Xmx{Xmx}g " \
                           "threads={threads} " \
                           "in={mpInterleavedDir}{sampleName}.{suffix} " \
                           "out={mp_split_read_folder}{sampleName}_P%_R1.fastq " \
                           "out2={mp_split_read_folder}{sampleName}_P%_R2.fastq " \
                           "ways={parts} 2>{split_log_folder}{sampleName}_split.log " \
                           .format(sampleName=self.sampleName,
                    suffix=GlobalParameter().suffix,
                    mp_split_read_folder=mp_split_read_folder,
                    split_log_folder=split_log_folder,
                    Xmx=GlobalParameter().maxMemory,
                    threads=GlobalParameter().threads,
                    parts=self.parts,
                    mpInterleavedDir=os.path.join(GlobalParameter().mpInterleavedDir + "/"))


        if self.libType == "pe":
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_pe)
            print(run_cmd(cmd_split_pe))

        if self.libType == "mp":
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_mp)
            print(run_cmd(cmd_split_mp))

        if self.libType == "pe-mp":
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_pe)
            print(run_cmd(cmd_split_pe))
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_mp)
            print(run_cmd(cmd_split_mp))

        if self.libType == "ilpe":
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_pe_il)
            print(run_cmd(cmd_split_pe_il))

        if self.libType == "ilmp":
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_mp_il)
            print(run_cmd(cmd_split_mp_il))

        if self.libType == "ilpe-ilmp":
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_pe_il)
            print(run_cmd(cmd_split_pe_il))
            print("****** NOW RUNNING COMMAND ******: " + cmd_split_mp_il)
            print(run_cmd(cmd_split_mp_il))



class splitSamples(luigi.Task):
    projectName = luigi.Parameter(default="ReadManipulation")
    parts = luigi.IntParameter(description="Number of parts")
    libType = luigi.ChoiceParameter(description="Choose From['pe: paired-end','mp: mate-pair','pe-mp: paired-end and mate-pair',"
                                                "'ilpe: interleaved paired-end','ilmp: interleaved mate-pair',"
                                                "'ilpe-ilmp: interleaved paired-end and interleaved mate-pair']",
                                    choices=["pe", "mp", "pe-mp", "ilpe", "ilmp", "ilpe-ilmp"], var_type=str)



    def requires(self):
        if self.libType == "pe":
            return [split(libType=self.libType,parts=self.parts,
                        sampleName=i)
                    for i in [line.strip()
                              for line in
                              open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
        if self.libType == "mp":
            return [split(libType=self.libType,parts=self.parts,
                           sampleName=i)
                    for i in [line.strip()
                              for line in
                              open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
        if self.libType == "pe-mp":
            return [
                        [split(libType="pe",parts=self.parts,sampleName=i)
                                for i in [line.strip()
                                          for line in
                                                open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                        [split(libType="mp",parts=self.parts, sampleName=i)
                                for i in [line.strip()
                                          for line in
                                                open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
                  ]

        if self.libType == "ilpe":
            return [split(libType=self.libType,parts=self.parts,
                        sampleName=i)
                    for i in [line.strip()
                              for line in
                              open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]]
        if self.libType == "ilmp":
            return [split(libType=self.libType,parts=self.parts,
                           sampleName=i)
                    for i in [line.strip()
                              for line in
                              open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]
        if self.libType == "ilpe-ilmp":
            return [
                [split(libType="ilpe",parts=self.parts, sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_il_samples.lst")))]],

                [split(libType="ilmp", parts=self.parts,sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_il_samples.lst")))]]
            ]


    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('task.splitsample.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('splitting finished at {t}'.format(t=timestamp))