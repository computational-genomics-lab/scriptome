import luigi
import os
import subprocess
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import filtlong

class GlobalParameter(luigi.Config):
    long_read_dir = luigi.Parameter()
    genome_size=luigi.Parameter()
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



def flye_formater(lrfile):
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_Long_Reads" + "/")

    with open(lrfile) as fh:
        sample_name_list = fh.read().splitlines()
        read_name_suffix = '.fastq'
        read_name_list = [lr_cleaned_read_folder + x + read_name_suffix for x in sample_name_list]
        lr_parse_string = ' '.join(read_name_list)
        return lr_parse_string


class flye(luigi.Task):
    projectName = luigi.Parameter(default="GenomeAssembly")

    read_library_type = luigi.ChoiceParameter(description="Choose From['pacbio raw, pacbio corrected, nanopore raw, "
                                                "nanopore corrected']",
                                    choices=["pacbio-raw", "pacbio-corr", "nano-raw", "nano-corr"], var_type=str)

    genome_size = luigi.Parameter(description="Estimated Genome Size in mb")

    threads=GlobalParameter().threads



    def requires(self):
        if self.read_library_type=="nano-raw" or self.read_library_type=="pacbio-raw":
            return [filtlong(sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]

    def output(self):
        flye_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "Flye_" + self.read_library_type +"/")
        return {'out': luigi.LocalTarget(flye_assembly_folder + "assembly.fasta")}

    def run(self):

        flye_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "Flye_" + self.read_library_type + "/")
        flye_assembly_log_folder = os.path.join(os.getcwd(), "log", "GenomeAssembly", "Flye_" + self.read_library_type +  "/")
        lr_sample_list = os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")
        cmd_flye_lr = flye_formater(lr_sample_list)



        flye_cmd = "[ -d  {flye_assembly_folder} ] || mkdir -p {flye_assembly_folder}; " \
                        "mkdir -p {flye_assembly_log_folder}; cd {flye_assembly_folder}; " \
                        "/usr/bin/time -v flye " \
                        "--{read_library_type} " \
                        "{cmd_flye_lr} " \
                        "--threads {threads} " \
                   "--genome-size {genomeSize}m " \
                   "--out-dir {flye_assembly_folder} " \
                   "2>&1 | tee {flye_assembly_log_folder}flye_assembly.log " \
            .format(flye_assembly_folder=flye_assembly_folder,
                    read_library_type=self.read_library_type,
                    cmd_flye_lr=cmd_flye_lr,
                    flye_assembly_log_folder=flye_assembly_log_folder,
                    genomeSize=self.genomeSize,
                    threads=GlobalParameter().threads)

        print("****** NOW RUNNING COMMAND ******: " + flye_cmd)
        run_cmd(flye_cmd)


