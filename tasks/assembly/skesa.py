#!/usr/bin/env python3
import os
from tasks.assembly.kmergenie import kmergenie_formater_reformat
from tasks.assembly.kmergenie import kmergenie_formater_bbduk
from tasks.assembly.kmergenie import optimal_kmer
from tasks.readCleaning.preProcessReads import bbduk
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


pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")


def skesa_pe_bbduk(samplefile):
    with open(samplefile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'
        left_read_name_list = [x + left_read_name_suffix for x in sample_name_list]
        right_read_name_list = [x + right_read_name_suffix for x in sample_name_list]
        pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
        result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]
        result1 = ["--fastq " + pe_cleaned_read_folder + x +"," + pe_cleaned_read_folder +y for x, y in result]
        parse_string = ' '.join(result1)
        return parse_string

def skesa_pe_reformat(samplefile):
    with open(samplefile) as fh:
        sample_name_list = fh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'
        left_read_name_list = [x + left_read_name_suffix for x in sample_name_list]
        right_read_name_list = [x + right_read_name_suffix for x in sample_name_list]
        pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
        result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]
        result1 = ["--fastq " + pe_cleaned_read_folder + x +"," + pe_cleaned_read_folder +y for x, y in result]
        parse_string = ' '.join(result1)
        return parse_string


class skesa(luigi.Task):
    projectName = luigi.Parameter(default="GenomeAssembly")
    domain=luigi.Parameter(default="prokaryote",description="Domian of the Organism. Must be Prokaryote Only")
    read_library_type = luigi.Parameter(default="pe")
    assembly_name = luigi.Parameter(default="assembly", description="Name of the Assembly")
    pre_process_reads = luigi.ChoiceParameter(choices=["yes","no"],var_type=str)
    kmer=luigi.IntParameter(default=21,description="Minimal Kmer Length for assembly. [--kmer 21]")
    steps=luigi.IntParameter(default=11,description="Number of assembly iterations from minimal to  maximal kmer length in reads [--steps 11]")
    min_contig_length=luigi.IntParameter(default=200,description="Minimal contig length reported in output [--min-contig-length 200]")

    def requires(self):

        if self.pre_process_reads=="yes":
            return [bbduk(read_library_type=self.read_library_type,
                      sampleName=i)
                for i in [line.strip()
                          for line in
                          open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]

        if self.pre_process_reads=="no":
            return [reformat(read_library_type=self.read_library_type,
                      sampleName=i)
                for i in [line.strip()
                          for line in
                          open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]

    def output(self):
        skesa_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "SKESA" + "/")
        return {'out': luigi.LocalTarget(skesa_assembly_folder + self.assembly_name + "-contigs.fa")}

    def run(self):
        skesa_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "SKESA" + "/")
        skesa_assembly_log_folder = os.path.join(os.getcwd(), "log", "GenomeAssembly", "SKESA" + "/")
        pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")

        if self.pre_process_reads=="yes":
            cmd_skesa_pe = skesa_pe_bbduk(pe_sample_list)
        if self.pre_process_reads=="no":
            cmd_skesa_pe = skesa_pe_reformat(pe_sample_list)

        run_cmd_skesa_pe = "[ -d  {skesa_assembly_folder} ] || mkdir -p {skesa_assembly_folder}; " \
                           "mkdir -p {skesa_assembly_log_folder}; cd {skesa_assembly_folder}; " \
                           "/usr/bin/time -v skesa " \
                           "--kmer {kmer} " \
                           "--steps {steps} --min_contig {min_contig_length} " \
                           "--cores {threads} " \
                           "--memory {maxMemory} " \
                           "--contigs_out {assembly_name}-contigs.fa " \
                           "{cmd_skesa_pe} " \
                           "2>&1 | tee {skesa_assembly_log_folder}skesa_assembly.log " \
            .format(skesa_assembly_folder=skesa_assembly_folder,
                    threads=GlobalParameter().threads,
                    kmer=self.kmer,
                    steps=self.steps,
                    min_contig_length=self.min_contig_length,
                    maxMemory=GlobalParameter().maxMemory,
                    assembly_name=self.assembly_name,
                    skesa_assembly_log_folder=skesa_assembly_log_folder,
                    cmd_skesa_pe=cmd_skesa_pe)
        if self.domain=="prokaryote":
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_skesa_pe)
            run_cmd(run_cmd_skesa_pe)