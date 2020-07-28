#!/usr/bin/env python3
import os
import re
from tasks.assembly.kmergenie import kmergenie_formater_bbduk
from tasks.assembly.kmergenie import kmergenie_formater_reformat
from tasks.assembly.kmergenie import optimal_kmer
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import filtlong
from tasks.readCleaning.reFormatReads import reformat


import luigi
import os
import subprocess

class GlobalParameter(luigi.Config):
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

minia=os.path.join(os.getcwd(),"GenomeAssembly", "MINIA")

createFolder(minia)


def minia_pe_bbduk(pefile):
    with open(pefile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'
        pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
        left_read_name_list = [ pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        right_read_name_list =[pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
        pe_reads = left_read_name_list + right_read_name_list

        with open((os.path.join(os.getcwd(),"GenomeAssembly", "MINIA","minia.fofn")), 'w') as fh:
            fh.writelines("%s\n"  % read for read in pe_reads)

def minia_pe_reformat(pefile):
    with open(pefile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'
        pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
        left_read_name_list = [ pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        right_read_name_list =[pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
        pe_reads = left_read_name_list + right_read_name_list

        with open((os.path.join(os.getcwd(),"GenomeAssembly", "MINIA","minia.fofn")), 'w') as fh:
            fh.writelines("%s\n"  % read for read in pe_reads)

def minia_pe_mp_bbduk(pefile,mpfile):

    with open(pefile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'

        pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")

        pe_left_read_name_list =  [ pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        pe_right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]

    with open(mpfile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'
        mp_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_MP_Reads" + "/")
        mp_left_read_name_list = [mp_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        mp_right_read_name_list = [mp_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]


        reads = pe_left_read_name_list + pe_right_read_name_list + mp_left_read_name_list + mp_right_read_name_list


        with open((os.path.join(os.getcwd(),"GenomeAssembly", "MINIA","minia.fofn")), 'w') as fh:
            fh.writelines("%s\n"  % read for read in reads)

def minia_pe_mp_reformat(pefile,mpfile):

    with open(pefile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'

        pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")

        pe_left_read_name_list =  [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        pe_right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]

    with open(mpfile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'
        mp_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_MP_Reads" + "/")
        mp_left_read_name_list = [mp_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        mp_right_read_name_list = [mp_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]


        reads = pe_left_read_name_list + pe_right_read_name_list + mp_left_read_name_list + mp_right_read_name_list


        with open((os.path.join(os.getcwd(),"GenomeAssembly", "MINIA","minia.fofn")), 'w') as fh:
            fh.writelines("%s\n"  % read for read in reads)
def minia_kmer(readlist):
    kmergenie_folder = os.path.join(os.getcwd(), "GenomeAssembly","MINIA" + "/")

    command = "[ -d  {kmergenie_folder} ] || mkdir -p {kmergenie_folder}; " \
                  "cd {kmergenie_folder}; " \
                  "kmergenie {kmergenie_folder}minia.fofn " \
                  .format(kmergenie_folder=kmergenie_folder)

    print("Estimating Optimal Kmers using kmergenie")
    print("-----------------------------------------")
    print("Command: ", command)
    print("\n")

    proc = subprocess.Popen(command,
                            shell=True,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    stdout_value, stderr_value = proc.communicate()

    parse_output = stdout_value.strip().decode("utf-8")

    p = re.compile(r'^best k:\s+(\d{2})', re.M)
    optimal_k = ' '.join(re.findall(p, parse_output))

    return optimal_k


class minia(luigi.Task):
    projectName = luigi.Parameter(default="GenomeAssembly")
    pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)

    read_library_type = luigi.ChoiceParameter(description="Choose From['pe-lr: paired-end and long read',"
                                                "'pe-mp-lr: paired-end, mate-pair and long read'",
                                    choices=["pe-lr","pe-mp-lr"], var_type=str)


    def requires(self):

        if self.read_library_type == "pe-lr" and self.pre_process_reads=="yes":
            return [
                [bbduk(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [minia_pe_bbduk(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [bbduk(read_library_type="lr",sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]

        if self.read_library_type == "pe-lr" and self.pre_process_reads=="no":
            return [
                [reformat(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [minia_pe_reformat(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [reformat(read_library_type="lr",sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]

        if self.read_library_type == "pe-mp-lr" and self.pre_process_reads=="yes":
            return [
                [bbduk(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],
                [bbduk(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]],

                [minia_pe_mp_bbduk((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")),
                                      (os.path.join(os.getcwd(), "sample_list", "mp_samples.lst"))
                                      )],

                [bbduk(read_library_type="lr",sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]

        if self.read_library_type == "pe-mp-lr" and self.pre_process_reads=="no":
            return [
                [reformat(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],
                [reformat(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]],

                [minia_pe_mp_reformat((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")),
                                      (os.path.join(os.getcwd(), "sample_list", "mp_samples.lst"))
                                      )],

                [bbduk(read_library_type="lr",sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]

    def output(self):
        minia_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "MINIA" + "/")
        return {'out': luigi.LocalTarget(minia_assembly_folder + "minia.contigs.fa")}

    def run(self):
        minia_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "MINIA" + "/")
        minia_assembly_log_folder = os.path.join(os.getcwd(), "log", "GenomeAssembly", "MINIA" +  "/")
        kmer = minia_kmer((os.path.join(os.getcwd(),"GenomeAssembly", "MINIA","minia.fofn")))
        print("Optimal Kmer: ", kmer)
        

        run_cmd_minia = "[ -d  {minia_assembly_folder} ] || mkdir -p {minia_assembly_folder}; " \
                        "mkdir -p {minia_assembly_log_folder}; cd {minia_assembly_folder}; " \
                        "/usr/bin/time -v minia " \
                        "-kmer-size {kmer} " \
                        "-nb-cores {threads} " \
                        "-in {minia_assembly_folder}minia.fofn " \
                        "2>&1 | tee {minia_assembly_log_folder}minia_assembly.log " \
            .format(minia_assembly_folder=minia_assembly_folder,
                    kmer=kmer,
                    threads=GlobalParameter().threads,
                    minia_assembly_log_folder=minia_assembly_log_folder)

        if self.read_library_type=="pe-lr" or self.read_library_type=="pe-mp-lr":

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_minia)
            print(run_cmd(run_cmd_minia))




