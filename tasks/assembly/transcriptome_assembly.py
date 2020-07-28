import luigi
import os
import time
import subprocess
import pandas as pd
from luigi import Parameter
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.reFormatReads import reformat
from tasks.assembly.spades import spades
from tasks.assembly.rockhopper import rockhopper
from tasks.assembly.denovo_trinity_assembly import trinity


def run_cmd(cmd):
    p = subprocess.Popen(cmd, bufsize=-1,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
    output = p.communicate()[0]
    return output


class GlobalParameter(luigi.Config):
    genome_suffix = luigi.Parameter()
    read_library_type = luigi.Parameter()
    organism_domain = luigi.Parameter()
    genome_name = luigi.Parameter()
    genome_dir = luigi.Parameter()
    transcriptome_dir = luigi.Parameter()
    threads = luigi.Parameter()
    maxMemory = luigi.Parameter()
    feature_type = luigi.Parameter()
    adapter = luigi.Parameter()


class DTA(luigi.Task):
    project_name = luigi.Parameter(default="RNASeqAnalysis")
    adapter = GlobalParameter().adapter
    organism_domain = GlobalParameter().organism_domain
    threads = GlobalParameter().threads
    pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
    read_library_type = GlobalParameter().read_library_type
    rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)

    def requires(self):
        if all([self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
            return [rockhopper(project_name=self.project_name,
                               pre_process_reads=self.pre_process_reads)
                    ]

        if self.rnaseq_assembler == "spades":
            return [spades(project_name="RNASeqAnalysis", read_library_type=self.read_library_type,
                           mode="rna",
                           pre_process_reads=self.pre_process_reads)
                    ]

        if self.rnaseq_assembler == "trinity":
            return [trinity(project_name=self.project_name,
                            pre_process_reads=self.pre_process_reads)
                    ]

    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget(os.path.join(os.getcwd(), "task_logs", 'task.assemble.transcript.complete.{t}'.format(t=timestamp)))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('transcript assembly finished at {t}'.format(t=timestamp))
