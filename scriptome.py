#!/usr/bin/env python3
import luigi
import sys

from tasks.dataDownload.fetch_data_ena import ENA
from tasks.dataDownload.fetch_data_ena import getDataENA
from tasks.dataDownload.fetch_data_ncbi import NCBI

from tasks.dataDownload.fetch_data_ncbi import getDataNCBI
from tasks.readManipulation.re_pair import rePair
from tasks.readManipulation.re_pair import rePairSamples
from tasks.readManipulation.verify_pairing import verifyPair
from tasks.readManipulation.verify_pairing import verifyPairSamples
from tasks.readManipulation.de_duplicate import deDup
from tasks.readManipulation.de_duplicate import deDupSamples
from tasks.readManipulation.fq_to_fasta import fq2Fasta
from tasks.readManipulation.fq_to_fasta import convertFq2Fasta
from tasks.readManipulation.il_to_pe import ilToPe
from tasks.readManipulation.il_to_pe import convertIlToPe
from tasks.readManipulation.sample import byNumber
from tasks.readManipulation.sample import sampleByNumber
from tasks.readManipulation.sample import byFraction
from tasks.readManipulation.sample import sampleByFraction
from tasks.readManipulation.split import split
from tasks.readManipulation.split import splitSamples

from tasks.readCleaning.rawReadQC import readqc
from tasks.readCleaning.rawReadQC import rawReadsQC
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import cleanReads
from tasks.readCleaning.preProcessReads import filtlong


from tasks.assembly import skesa
from tasks.assembly import redundans
from tasks.assembly import flye
from tasks.assembly import minia
from tasks.assembly import dbg2olc
from tasks.assembly import rockhopper
from tasks.assembly.spades import *
from tasks.assembly.unicycler import *
from tasks.assembly.kmergenie import *
from tasks.assembly import abyss
#from tasks.assembly import metam_16_july
from tasks.assembly import metagenome_analysis
from tasks.assembly import meta_profile
from tasks.assembly import denovo_trinity_assembly
from tasks.assembly.transcriptome_assembly import DTA


from tasks.rnaSeq import index_genome
from tasks.rnaSeq import align_rnaseq_reads_with_genome
from tasks.rnaSeq import index_transctriptome
from tasks.rnaSeq import index_quanify_cluster_da_transcripts
from tasks.rnaSeq import index_quanify_cluster_gg_transcripts
from tasks.rnaSeq import generate_gene_count_file
from tasks.rnaSeq import generate_transcript_count_file
from tasks.rnaSeq import denovo_differential_expression_analysis
from tasks.rnaSeq import genome_guided_differential_expression_analysis
from tasks.rnaSeq import alignment_based_differential_expression_analysis
from tasks.rnaSeq import alignment_free_differential_expression_analysis


from tasks.annotation import prokaryotic_annotation
from tasks.annotation.make_tx_to_gene import makeTx2Gene


if __name__ == '__main__':

    luigi.run()
