#!/usr/bin/env python3
import os
import re
from tasks.assembly.kmergenie import kmergenie_formater_bbduk
from tasks.assembly.kmergenie import kmergenie_formater_reformat
from tasks.assembly.kmergenie import optimal_kmer
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import filtlong
from tasks.assembly.minia import *

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

#minia=os.path.join(os.getcwd(),"GenomeAssembly", "MINIA")

#createFolder(minia)




def dbg2olc_formater(lrfile):
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_Long_Reads" + "/")

    with open(lrfile) as fh:
        sample_name_list = fh.read().splitlines()
        read_name_suffix = '.fastq'
        read_name_list = ["f " + lr_cleaned_read_folder + x + read_name_suffix for x in sample_name_list]
        lr_parse_string =  ' '.join(read_name_list)
        return lr_parse_string



class dbg2olc(luigi.Task):
    projectName = luigi.Parameter(default="GenomeAssembly")
    pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
    read_library_type = luigi.ChoiceParameter(description="Choose From['pe-lr: paired-end and long read',"
                                                "'pe-mp-lr: paired-end, mate-pair and long read'",
                                    choices=["pe-lr","pe-mp-lr"], var_type=str)


    def requires(self):

        if self.read_library_type == "pe-lr" or self.read_library_type == "pe-mp-lr":
            return [minia(read_library_type=self.read_library_type,
                          pre_process_reads=self.pre_process_reads)]


    def output(self):
        dbg2olc_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "DBG2OLC" + "/")
        return {'out': luigi.LocalTarget(dbg2olc_assembly_folder + "DBG2OLC_contigs.fa")}

    def run(self):
        minia_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "MINIA" + "/")
        dbg2olc_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "DBG2OLC" + "/")

        DBG2OLC_assembly_log_folder = os.path.join(os.getcwd(), "log", "GenomeAssembly", "DBG2OLC" +  "/")

        kmer = minia_kmer((os.path.join(os.getcwd(),"GenomeAssembly", "MINIA","minia.fofn")))

        print("Optimal Kmer: ", kmer)

        dbg2olc_input=dbg2olc_formater((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))

        run_cmd_dbg2olc = "[ -d  {dbg2olc_assembly_folder} ] || mkdir -p {dbg2olc_assembly_folder}; " \
                        "mkdir -p {DBG2OLC_assembly_log_folder}; cd {dbg2olc_assembly_folder}; " \
                        "/usr/bin/time -v DBG2OLC " \
                        "k {kmer} Contigs {minia_assembly_folder}minia.contigs.fa " \
                        "KmerCovTh 2 MinOverlap 20 AdaptiveTh 0.005 " \
                        "{dbg2olc_input} " \
                        "2>&1 | tee {DBG2OLC_assembly_log_folder}dbg2olc_assembly.log " \
            .format(minia_assembly_folder=minia_assembly_folder,
                    dbg2olc_assembly_folder=dbg2olc_assembly_folder,
                    dbg2olc_input=dbg2olc_input,
                    DBG2OLC_assembly_log_folder=DBG2OLC_assembly_log_folder,
                    kmer=kmer)

        if self.read_library_type == "pe-lr" or self.read_library_type == "pe-mp-lr":

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_dbg2olc)
            print(run_cmd(run_cmd_dbg2olc))




