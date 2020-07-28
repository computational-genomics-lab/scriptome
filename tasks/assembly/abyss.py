#!/usr/bin/env python3
import os
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


pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")

def abyss_pe_bbduk(samplefile):
    with open(samplefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]

        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'"

        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'

        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]

        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]

        pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")

        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

        pe_result2 = ' '.join(pe_result1)

        pe_parse_string = pe_lib_string + pe_result2
        return pe_parse_string


def abyss_pe_mp_bbduk(pefile, mpfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_MP_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)

        parse_string = pe_lib_string + mp_lib_string + pe_result2 + mp_result2

        return parse_string


def abyss_pe_mp_se_bbduk(pefile, mpfile, sefile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_MP_Reads" + "/")
    se_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_SE_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

    with open(sefile) as fh:
        se_sample_name_list = fh.read().splitlines()
        se_samples = [x for x in se_sample_name_list]
        se_lib_string = "se=" + "'" + ' '.join(se_samples) + "'" + " "
        se_read_name_suffix = '.fastq'
        se_read_name_list = [x + se_read_name_suffix for x in se_sample_name_list]
        se_result = [sublist for sublist in zip(se_samples, se_read_name_list)]
        se_result1 = [' ' + x + '=' + "'" + se_cleaned_read_folder + y + "'" for x, y in se_result]

        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)
        se_result2 = ' '.join(se_result1)

        parse_string = pe_lib_string + mp_lib_string + se_lib_string + pe_result2 + mp_result2 + se_result2

        return parse_string


def abyss_pe_mp_se_lr_bbduk(pefile, mpfile, sefile, lrfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_MP_Reads" + "/")
    se_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_SE_Reads" + "/")
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_Long_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

    with open(sefile) as fh:
        se_sample_name_list = fh.read().splitlines()
        se_samples = [x for x in se_sample_name_list]
        se_lib_string = "se=" + "'" + ' '.join(se_samples) + "'" + " "
        se_read_name_suffix = '.fastq'
        se_read_name_list = [x + se_read_name_suffix for x in se_sample_name_list]
        se_result = [sublist for sublist in zip(se_samples, se_read_name_list)]
        se_result1 = [' ' + x + '=' + "'" + se_cleaned_read_folder + y + "'" for x, y in se_result]

    with open(lrfile) as fh:
        lr_sample_name_list = fh.read().splitlines()
        lr_samples = [x for x in lr_sample_name_list]
        lr_lib_string = "long=" + "'" + ' '.join(lr_samples) + "'" + " "
        lr_read_name_suffix = '.fastq'
        lr_read_name_list = [x + lr_read_name_suffix for x in lr_sample_name_list]
        lr_result = [sublist for sublist in zip(lr_samples, lr_read_name_list)]
        lr_result1 = [' ' + x + '=' + "'" + lr_cleaned_read_folder + y + "'" for x, y in lr_result]

        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)
        se_result2 = ' '.join(se_result1)
        lr_result2 = ' '.join(lr_result1)

        parse_string = pe_lib_string + mp_lib_string + se_lib_string + lr_lib_string + pe_result2 + mp_result2 + se_result2 + lr_result2

        return parse_string


def abyss_pe_mp_lr_bbduk(pefile, mpfile, lrfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_MP_Reads" + "/")
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_Long_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

    with open(lrfile) as fh:
        lr_sample_name_list = fh.read().splitlines()
        lr_samples = [x for x in lr_sample_name_list]
        lr_lib_string = "long=" + "'" + ' '.join(lr_samples) + "'" + " "
        lr_read_name_suffix = '.fastq'
        lr_read_name_list = [x + lr_read_name_suffix for x in lr_sample_name_list]
        lr_result = [sublist for sublist in zip(lr_samples, lr_read_name_list)]
        lr_result1 = [' ' + x + '=' + "'" + lr_cleaned_read_folder + y + "'" for x, y in lr_result]

        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)
        lr_result2 = ' '.join(lr_result1)

        parse_string = pe_lib_string + mp_lib_string + lr_lib_string + pe_result2 + mp_result2 + lr_result2

        return parse_string


###########PE with LR#######################
############################################
def abyss_pe_lr_bbduk(pefile, lrfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_Long_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(lrfile) as fh:
        lr_sample_name_list = fh.read().splitlines()
        lr_samples = [x for x in lr_sample_name_list]
        lr_lib_string = "long=" + "'" + ' '.join(lr_samples) + "'" + " "
        lr_read_name_suffix = '.fastq'
        lr_read_name_list = [x + lr_read_name_suffix for x in lr_sample_name_list]
        lr_result = [sublist for sublist in zip(lr_samples, lr_read_name_list)]
        lr_result1 = [' ' + x + '=' + "'" + lr_cleaned_read_folder + y + "'" for x, y in lr_result]

        pe_result2 = ' '.join(pe_result1)
        lr_result2 = ' '.join(lr_result1)

        parse_string = pe_lib_string + lr_lib_string + pe_result2 + lr_result2

        return parse_string
#####################################################################################################################################
#RUN REFORMAT.SH IF NO READ CLEANING
#####################################################################################################################################
def abyss_pe_reformat(samplefile):
    with open(samplefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]

        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'"

        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'

        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]

        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]

        pe_cleaned_read_folder = os.path.join(os.getcwd(),  "VerifiedReads", "Verified_PE_Reads" + "/")

        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

        pe_result2 = ' '.join(pe_result1)

        pe_parse_string = pe_lib_string + pe_result2
        return pe_parse_string


def abyss_pe_mp_reformat(pefile, mpfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_MP_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)

        parse_string = pe_lib_string + mp_lib_string + pe_result2 + mp_result2

        return parse_string

def abyss_pe_mp_lr_reformat(pefile, mpfile, lrfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_MP_Reads" + "/")
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_Long_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

    with open(lrfile) as fh:
        lr_sample_name_list = fh.read().splitlines()
        lr_samples = [x for x in lr_sample_name_list]
        lr_lib_string = "long=" + "'" + ' '.join(lr_samples) + "'" + " "
        lr_read_name_suffix = '.fastq'
        lr_read_name_list = [x + lr_read_name_suffix for x in lr_sample_name_list]
        lr_result = [sublist for sublist in zip(lr_samples, lr_read_name_list)]
        lr_result1 = [' ' + x + '=' + "'" + lr_cleaned_read_folder + y + "'" for x, y in lr_result]

        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)
        lr_result2 = ' '.join(lr_result1)

        parse_string = pe_lib_string + mp_lib_string + lr_lib_string + pe_result2 + mp_result2 + lr_result2

        return parse_string


def abyss_pe_mp_se_reformat(pefile, mpfile, sefile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_MP_Reads" + "/")
    se_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_SE_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

    with open(sefile) as fh:
        se_sample_name_list = fh.read().splitlines()
        se_samples = [x for x in se_sample_name_list]
        se_lib_string = "se=" + "'" + ' '.join(se_samples) + "'" + " "
        se_read_name_suffix = '.fastq'
        se_read_name_list = [x + se_read_name_suffix for x in se_sample_name_list]
        se_result = [sublist for sublist in zip(se_samples, se_read_name_list)]
        se_result1 = [' ' + x + '=' + "'" + se_cleaned_read_folder + y + "'" for x, y in se_result]
        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)
        se_result2 = ' '.join(se_result1)
        parse_string = pe_lib_string + mp_lib_string + se_lib_string + pe_result2 + mp_result2 + se_result2
        return parse_string


def abyss_pe_mp_se_lr_reformat(pefile, mpfile, sefile, lrfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
    mp_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_MP_Reads" + "/")
    se_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_SE_Reads" + "/")
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_Long_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(mpfile) as fh:
        mp_sample_name_list = fh.read().splitlines()
        mp_samples = [x for x in mp_sample_name_list]
        mp_lib_string = "mp=" + "'" + ' '.join(mp_samples) + "'" + " "
        mp_left_read_name_suffix = '_R1.fastq'
        mp_right_read_name_suffix = '_R2.fastq'
        mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
        mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
        mp_samples = [x for x in mp_sample_name_list]
        mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
        mp_result1 = [' ' + x + '=' + "'" + mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for
                      x, y, z in mp_result]

    with open(sefile) as fh:
        se_sample_name_list = fh.read().splitlines()
        se_samples = [x for x in se_sample_name_list]
        se_lib_string = "se=" + "'" + ' '.join(se_samples) + "'" + " "
        se_read_name_suffix = '.fastq'
        se_read_name_list = [x + se_read_name_suffix for x in se_sample_name_list]
        se_result = [sublist for sublist in zip(se_samples, se_read_name_list)]
        se_result1 = [' ' + x + '=' + "'" + se_cleaned_read_folder + y + "'" for x, y in se_result]

    with open(lrfile) as fh:
        lr_sample_name_list = fh.read().splitlines()
        lr_samples = [x for x in lr_sample_name_list]
        lr_lib_string = "long=" + "'" + ' '.join(lr_samples) + "'" + " "
        lr_read_name_suffix = '.fastq'
        lr_read_name_list = [x + lr_read_name_suffix for x in lr_sample_name_list]
        lr_result = [sublist for sublist in zip(lr_samples, lr_read_name_list)]
        lr_result1 = [' ' + x + '=' + "'" + lr_cleaned_read_folder + y + "'" for x, y in lr_result]

        pe_result2 = ' '.join(pe_result1)
        mp_result2 = ' '.join(mp_result1)
        se_result2 = ' '.join(se_result1)
        lr_result2 = ' '.join(lr_result1)

        parse_string = pe_lib_string + mp_lib_string + se_lib_string + lr_lib_string + pe_result2 + mp_result2 + se_result2 + lr_result2

        return parse_string


###########PE with LR#######################
############################################
def abyss_pe_lr_reformat(pefile, lrfile):
    pe_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
    lr_cleaned_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_Long_Reads" + "/")

    with open(pefile) as fh:
        pe_sample_name_list = fh.read().splitlines()
        pe_samples = [x for x in pe_sample_name_list]
        pe_lib_string = "lib=" + "'" + ' '.join(pe_samples) + "'" + " "
        pe_left_read_name_suffix = '_R1.fastq'
        pe_right_read_name_suffix = '_R2.fastq'
        pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
        pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
        pe_samples = [x for x in pe_sample_name_list]
        pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
        pe_result1 = [' ' + x + '=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for
                      x, y, z in pe_result]

    with open(lrfile) as fh:
        lr_sample_name_list = fh.read().splitlines()
        lr_samples = [x for x in lr_sample_name_list]
        lr_lib_string = "long=" + "'" + ' '.join(lr_samples) + "'" + " "
        lr_read_name_suffix = '.fastq'
        lr_read_name_list = [x + lr_read_name_suffix for x in lr_sample_name_list]
        lr_result = [sublist for sublist in zip(lr_samples, lr_read_name_list)]
        lr_result1 = [' ' + x + '=' + "'" + lr_cleaned_read_folder + y + "'" for x, y in lr_result]

        pe_result2 = ' '.join(pe_result1)
        lr_result2 = ' '.join(lr_result1)

        parse_string = pe_lib_string + lr_lib_string + pe_result2 + lr_result2

        return parse_string
#####################################################################################################################################
class abyss(luigi.Task):
    projectName = luigi.Parameter(default="GenomeAssembly")
    pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
    read_library_type = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair'",
                                              choices=["pe", "pe-mp", "pe-lr", "pe-mp-lr", "pe-mp-se", "pe-mp-se-lr"], var_type=str)

    assembly_name = luigi.Parameter(default="assembly", description="Name of the Assembly")

    def requires(self):
        if self.read_library_type == "pe" and self.pre_process_reads =="yes":
            return [[bbduk(read_library_type=self.read_library_type,
                           sampleName=i)
                     for i in [line.strip()
                               for line in
                               open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                    [kmergenie_formater_bbduk(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],
                    ]

        if self.read_library_type == "pe" and self.pre_process_reads =="no":
            return [[reformat(read_library_type=self.read_library_type,
                           sampleName=i)
                     for i in [line.strip()
                               for line in
                               open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                    [kmergenie_formater_reformat(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],
                    ]

        if self.read_library_type == "pe-mp" and self.pre_process_reads =="yes":
            return [
                [bbduk(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_bbduk(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [bbduk(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
            ]

        if self.read_library_type == "pe-mp" and self.pre_process_reads =="no":
            return [
                [reformat(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_reformat(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [reformat(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
            ]




        if self.read_library_type == "pe-lr" and self.pre_process_reads == "yes":
            return [
                [bbduk(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_bbduk(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [filtlong(sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]

        if self.read_library_type == "pe-lr" and self.pre_process_reads == "no":
            return [
                [reformat(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_reformat(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [reformat(sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]

        if self.read_library_type == "pe-mp-se" and self.pre_process_reads == "yes":
            return [
                [bbduk(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_bbduk(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [bbduk(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]],

                [bbduk(read_library_type="se", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]
            ]

        if self.read_library_type == "pe-mp-se" and self.pre_process_reads == "no":
            return [
                [reformat(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_reformat(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [reformat(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]],

                [reformat(read_library_type="se", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "se_samples.lst")))]]
            ]

        if self.read_library_type == "pe-mp-lr" and self.pre_process_reads == "yes":
            return [
                [bbduk(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_bbduk(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [bbduk(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]],

                [bbduk(read_library_type="lr", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]
        if self.read_library_type == "pe-mp-lr" and self.pre_process_reads == "no":
            return [
                [reformat(read_library_type="pe", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                [kmergenie_formater_reformat(os.path.join(os.getcwd(), "sample_list", "pe_samples.lst"))],

                [reformat(read_library_type="mp", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]],

                [reformat(read_library_type="lr", sampleName=i)
                 for i in [line.strip()
                           for line in
                           open((os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")))]]
            ]

    def output(self):
        abyss_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "ABYSS" + "/")
        return {'out': luigi.LocalTarget(abyss_assembly_folder + self.assembly_name + "-scaffolds.fa")}

    def run(self):
        abyss_assembly_folder = os.path.join(os.getcwd(), "GenomeAssembly", "ABYSS" + "/")

        abyss_assembly_log_folder = os.path.join(os.getcwd(), "log", "GenomeAssembly", "ABYSS" + "/")

        pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")
        mp_sample_list = os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")
        se_sample_list = os.path.join(os.getcwd(), "sample_list", "se_samples.lst")
        lr_sample_list = os.path.join(os.getcwd(), "sample_list", "lr_samples.lst")


        if self.pre_process_reads == "yes":
            cmd_abyss_pe = abyss_pe_bbduk(pe_sample_list)
            cmd_abyss_pe_mp = abyss_pe_mp_bbduk(pe_sample_list,mp_sample_list)
            cmd_abyss_pe_mp_se = abyss_pe_mp_se_bbduk(pe_sample_list, mp_sample_list, se_sample_list)
            cmd_abyss_pe_mp_lr = abyss_pe_mp_se_bbduk(pe_sample_list, mp_sample_list, lr_sample_list)
            cmd_abyss_pe_lr = abyss_pe_lr_bbduk(pe_sample_list, lr_sample_list)
            cmd_abyss_pe_mp_se_lr = abyss_pe_mp_se_lr_bbduk(pe_sample_list, mp_sample_list, se_sample_list, lr_sample_list)

        if self.pre_process_reads == "no":
            cmd_abyss_pe = abyss_pe_reformat(pe_sample_list)
            cmd_abyss_pe_mp=abyss_pe_mp_reformat(pe_sample_list,mp_sample_list)
            cmd_abyss_pe_mp_lr = abyss_pe_mp_lr_reformat(pe_sample_list, mp_sample_list, lr_sample_list)
            cmd_abyss_pe_lr = abyss_pe_lr_reformat(pe_sample_list, lr_sample_list)
            cmd_abyss_pe_mp_se = abyss_pe_mp_se_reformat(pe_sample_list, mp_sample_list, se_sample_list)
            cmd_abyss_pe_mp_se_lr = abyss_pe_mp_se_lr_reformat(pe_sample_list, mp_sample_list, se_sample_list, lr_sample_list)



        kmer = optimal_kmer((os.path.join(os.getcwd(), "GenomeAssembly", "KmerGenie", "kmergenni_pe.lst")))

        run_cmd_abyss_pe = "[ -d  {abyss_assembly_folder} ] || mkdir -p {abyss_assembly_folder}; " \
                           "mkdir -p {abyss_assembly_log_folder}; cd {abyss_assembly_folder}; " \
                           "/usr/bin/time -v abyss-pe " \
                           "k={kmer} " \
                           "np={threads} " \
                           "name={assembly_name} " \
                           "{cmd_abyss_pe} " \
                           "-o {abyss_assembly_folder} " \
                           "2>&1 | tee {abyss_assembly_log_folder}abyss_assembly.log " \
            .format(abyss_assembly_folder=abyss_assembly_folder,
                    kmer=kmer,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name,
                    abyss_assembly_log_folder=abyss_assembly_log_folder,
                    cmd_abyss_pe=cmd_abyss_pe)

        run_cmd_abyss_pe_mp = "[ -d  {abyss_assembly_folder} ] || mkdir -p {abyss_assembly_folder}; " \
                              "mkdir -p {abyss_assembly_log_folder}; cd {abyss_assembly_folder}; " \
                              "/usr/bin/time -v abyss-pe " \
                              "k={kmer} " \
                              "np={threads} " \
                              "name={assembly_name} " \
                              "{cmd_abyss_pe_mp} " \
                              "-o {abyss_assembly_folder} " \
                              "2>&1 | tee {abyss_assembly_log_folder}abyss_assembly.log " \
            .format(abyss_assembly_folder=abyss_assembly_folder,
                    kmer=kmer,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name,
                    abyss_assembly_log_folder=abyss_assembly_log_folder,
                    cmd_abyss_pe_mp=cmd_abyss_pe_mp)

        run_cmd_abyss_pe_mp_lr = "[ -d  {abyss_assembly_folder} ] || mkdir -p {abyss_assembly_folder}; " \
                              "mkdir -p {abyss_assembly_log_folder}; cd {abyss_assembly_folder}; " \
                              "/usr/bin/time -v abyss-pe " \
                              "k={kmer} " \
                              "np={threads} " \
                              "name={assembly_name} " \
                              "{cmd_abyss_pe_mp_lr} " \
                              "-o {abyss_assembly_folder} " \
                              "2>&1 | tee {abyss_assembly_log_folder}abyss_assembly.log " \
            .format(abyss_assembly_folder=abyss_assembly_folder,
                    kmer=kmer,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name,
                    abyss_assembly_log_folder=abyss_assembly_log_folder,
                    cmd_abyss_pe_mp_lr=cmd_abyss_pe_mp_lr)

        run_cmd_abyss_pe_lr = "[ -d  {abyss_assembly_folder} ] || mkdir -p {abyss_assembly_folder}; " \
                              "mkdir -p {abyss_assembly_log_folder}; cd {abyss_assembly_folder}; " \
                              "/usr/bin/time -v abyss-pe " \
                              "k={kmer} " \
                              "np={threads} " \
                              "name={assembly_name} " \
                              "{cmd_abyss_pe_lr} " \
                              "-o {abyss_assembly_folder} " \
                              "2>&1 | tee {abyss_assembly_log_folder}abyss_assembly.log " \
            .format(abyss_assembly_folder=abyss_assembly_folder,
                    kmer=kmer,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name,
                    abyss_assembly_log_folder=abyss_assembly_log_folder,
                    cmd_abyss_pe_lr=cmd_abyss_pe_lr)

        run_cmd_abyss_pe_mp_se = "[ -d  {abyss_assembly_folder} ] || mkdir -p {abyss_assembly_folder}; " \
                                 "mkdir -p {abyss_assembly_log_folder}; cd {abyss_assembly_folder}; " \
                                 "/usr/bin/time -v abyss-pe " \
                                 "k={kmer} " \
                                 "np={threads} " \
                                 "name={assembly_name} " \
                                 "{cmd_abyss_pe_mp_se} " \
                                 "-o {abyss_assembly_folder} " \
                                 "2>&1 | tee {abyss_assembly_log_folder}abyss_assembly.log " \
            .format(abyss_assembly_folder=abyss_assembly_folder,
                    kmer=kmer,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name,
                    abyss_assembly_log_folder=abyss_assembly_log_folder,
                    cmd_abyss_pe_mp_se=cmd_abyss_pe_mp_se)

        run_cmd_abyss_pe_mp_se_lr = "[ -d  {abyss_assembly_folder} ] || mkdir -p {abyss_assembly_folder}; " \
                                 "mkdir -p {abyss_assembly_log_folder}; cd {abyss_assembly_folder}; " \
                                 "/usr/bin/time -v abyss-pe " \
                                 "k={kmer} " \
                                 "np={threads} " \
                                 "name={assembly_name} " \
                                 "{cmd_abyss_pe_mp_se_lr} " \
                                 "-o {abyss_assembly_folder} " \
                                 "2>&1 | tee {abyss_assembly_log_folder}abyss_assembly.log " \
            .format(abyss_assembly_folder=abyss_assembly_folder,
                    kmer=kmer,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name,
                    abyss_assembly_log_folder=abyss_assembly_log_folder,
                    cmd_abyss_pe_mp_se_lr=cmd_abyss_pe_mp_se_lr)

        if self.read_library_type == "pe" and self.pre_process_reads == "yes":
            kmergenie_formater_bbduk(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe)
            run_cmd(run_cmd_abyss_pe)

        if self.read_library_type == "pe" and self.pre_process_reads == "no":
            kmergenie_formater_reformat(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe)
            run_cmd(run_cmd_abyss_pe)

        if self.read_library_type == "pe-mp" and self.pre_process_reads == "yes":
            kmergenie_formater_bbduk(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp)
            run_cmd(run_cmd_abyss_pe_mp)

        if self.read_library_type == "pe-mp" and self.pre_process_reads == "no":
            kmergenie_formater_reformat(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp)
            run_cmd(run_cmd_abyss_pe_mp)

        if self.read_library_type == "pe-lr" and self.pre_process_reads == "yes":
            kmergenie_formater_bbduk(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_lr)
            run_cmd(run_cmd_abyss_pe_lr)

        if self.read_library_type == "pe-lr" and self.pre_process_reads == "no":
            kmergenie_formater_reformat(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_lr)
            run_cmd(run_cmd_abyss_pe_lr)

        if self.read_library_type == "pe-mp-se" and self.pre_process_reads=="yes":
            kmergenie_formater_bbduk(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp_se)
            run_cmd(run_cmd_abyss_pe_mp_se)

        if self.read_library_type == "pe-mp-se" and self.pre_process_reads=="no":
            kmergenie_formater_reformat(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp_se)
            run_cmd(run_cmd_abyss_pe_mp_se)

        if self.read_library_type == "pe-mp-lr" and self.pre_process_reads=="yes":
            kmergenie_formater_bbduk(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp_lr)
            run_cmd(run_cmd_abyss_pe_mp_lr)

        if self.read_library_type == "pe-mp-lr" and self.pre_process_reads=="no":
            kmergenie_formater_reformat(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp_lr)
            run_cmd(run_cmd_abyss_pe_mp_lr)

        if self.read_library_type == "pe-mp-se-lr" and self.pre_process_reads=="yes":
            kmergenie_formater_bbduk(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp_se_lr)
            run_cmd(run_cmd_abyss_pe_mp_se_lr)

        if self.read_library_type == "pe-mp-se-lr" and self.pre_process_reads=="yes":
            kmergenie_formater_bbduk(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp_se_lr)
            run_cmd(run_cmd_abyss_pe_mp_se_lr)
        if self.read_library_type == "pe-mp-se-lr" and self.pre_process_reads == "no":
            kmergenie_formater_reformat(pe_sample_list)
            print("Optimal Kmer: ", kmer)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_abyss_pe_mp_se_lr)
            run_cmd(run_cmd_abyss_pe_mp_se_lr)



