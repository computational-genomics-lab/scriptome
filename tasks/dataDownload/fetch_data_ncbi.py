import luigi
import os
import time
import subprocess

def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output


class NCBI(luigi.Task):
    #projectName = luigi.Parameter()

    section = luigi.Parameter(description="Choose from:['refseq','genbank']")

    fetch = luigi.Parameter(description="Choose from:['assembly-fna','assembly-gbk', 'assembly-report', 'assembly-stats', 'assembly-gff3', "
                                        "'protseq', 'nuclseq']")
    accession = luigi.Parameter("download sequences matching the provided NCBI-assembly accession(s) or NCBI Sequence ID"
                                                "A comma-separated list of accessions is possible, "
                                                "as well as a path to a filename containing one accession per line")
    def output(self):
        downloadFolder = os.path.join(os.getcwd(), "raw_data/")

        if self.fetch == "assembly-fna":
            return {'out': luigi.LocalTarget(downloadFolder  + self.accession + "_genomic.fna.gz")}
        if self.fetch == "assembly-gff3":
            return {'out': luigi.LocalTarget(downloadFolder  + self.accession + "_genomic.gff.gz")}


        if self.fetch == "assembly-gbk":
            return {'out': luigi.LocalTarget(downloadFolder  + self.accession + "_genomic.fna.gz")}
     
        if self.fetch == "assembly-report":
            return {'out': luigi.LocalTarget(downloadFolder  + self.accession + "_assembly_report.txt")}

        if self.fetch == "assembly-stats":
            return {'out': luigi.LocalTarget(downloadFolder  + self.accession + "_assembly_stats.txt")}


        if self.fetch == "protseq":
            return {'out': luigi.LocalTarget(downloadFolder  + self.accession+".fa")}
        if self.fetch == "nuclseq":
            return {'out': luigi.LocalTarget(downloadFolder  + self.accession+".fa")}

    def run(self):
        downloadFolder = os.path.join(os.getcwd(), "raw_data/")
        tmpFolder = os.path.join(os.getcwd(), "tmp")

        cmd_purge_tmp = "rm -rf {tmpFolder}".format(tmpFolder=tmpFolder)

        cmd_get_gff = "[ -d  {downloadFolder} ] || mkdir -p {downloadFolder}; " \
                      "[ -d  {tmpFolder} ] || mkdir -p {tmpFolder}; " \
                      "ncbi-genome-download -o {tmpFolder} -F gff -s {section} -A {accession} all" \
                .format(downloadFolder=downloadFolder, tmpFolder=tmpFolder, accession=self.accession, section=self.section)

        cp_gff = "find ./tmp -name \*gff.gz -exec cp {} "
        cmd_cp_gff = "{cp_gff} {downloadFolder}{accession}_genomic.gff.gz \;".format(cp_gff=cp_gff,
                                                                                     accession=self.accession,
                                                                                     downloadFolder=downloadFolder)


        cmd_get_assembly = "[ -d  {downloadFolder} ] || mkdir -p {downloadFolder}; " \
                           "[ -d  {tmpFolder} ] || mkdir -p {tmpFolder}; " \
                           "ncbi-genome-download -o {tmpFolder} -F fasta -s {section} -A {accession} all" \
                .format(downloadFolder=downloadFolder, tmpFolder=tmpFolder, accession=self.accession, section=self.section)

        cp_assembly = "find ./tmp -name \*fna.gz -exec cp {} "
        cmd_cp_assembly = "{cp_assembly} {downloadFolder}{accession}_genomic.fna.gz \;".format(cp_assembly=cp_assembly,
                                                                     downloadFolder=downloadFolder,accession=self.accession)

        cmd_get_assembly_gbk = "[ -d  {downloadFolder} ] || mkdir -p {downloadFolder}; " \
                           "[ -d  {tmpFolder} ] || mkdir -p {tmpFolder}; " \
                           "ncbi-genome-download -o {tmpFolder} -F genbank -s {section} -A {accession} all" \
                .format(downloadFolder=downloadFolder, tmpFolder=tmpFolder, accession=self.accession, section=self.section)

        cp_assembly_gbk = "find ./tmp -name \*gbff.gz -exec cp {} "
        cmd_cp_assembly_gbk = "{cp_assembly} {downloadFolder}{accession}_genomic.gbff.gz \;".format(cp_assembly=cp_assembly_gbk,
                                                                     downloadFolder=downloadFolder,accession=self.accession)

        cmd_get_assembly_report = "[ -d  {downloadFolder} ] || mkdir -p {downloadFolder}; " \
                           "[ -d  {tmpFolder} ] || mkdir -p {tmpFolder}; " \
                           "ncbi-genome-download -o {tmpFolder} -F assembly-report -s {section} -A {accession} all" \
                .format(downloadFolder=downloadFolder, tmpFolder=tmpFolder, accession=self.accession, section=self.section)

        cp_assembly_report = "find ./tmp -name \*assembly_report.txt -exec cp {} "
        cmd_cp_assembly_report = "{cp_assembly} {downloadFolder}{accession}_assembly_report.txt \;".format(cp_assembly=cp_assembly_report,
                                                                     downloadFolder=downloadFolder,accession=self.accession)

        cmd_get_assembly_stats = "[ -d  {downloadFolder} ] || mkdir -p {downloadFolder}; " \
                           "[ -d  {tmpFolder} ] || mkdir -p {tmpFolder}; " \
                           "ncbi-genome-download -o {tmpFolder} -F assembly-stats -s {section} -A {accession} all" \
                .format(downloadFolder=downloadFolder, tmpFolder=tmpFolder, accession=self.accession, section=self.section)

        cp_assembly_stats = "find ./tmp -name \*assembly_stats.txt -exec cp {} "
        cmd_cp_assembly_stats = "{cp_assembly} {downloadFolder}{accession}_assembly_stats.txt \;".format(cp_assembly=cp_assembly_report,
                                                                     downloadFolder=downloadFolder,accession=self.accession)




        cmd_get_nucleotide_fasta = "[ -d  {downloadFolder} ] || mkdir -p {downloadFolder}; cd {downloadFolder}; " \
                                 "ncbi-acc-download -m nucleotide -F fasta {accession} " \
                .format(downloadFolder=downloadFolder, accession=self.accession)

        cmd_get_protein_fasta = "[ -d  {downloadFolder} ] || mkdir -p {downloadFolder}; cd {downloadFolder}; " \
                                   "ncbi-acc-download -m protein -F fasta {accession} " \
            .format(downloadFolder=downloadFolder, accession=self.accession)


        if self.fetch == "assembly-gff3":
            print("****** NOW RUNNING COMMAND ******: " + cmd_get_gff)
            run_cmd(cmd_get_gff)
            print("****** NOW RUNNING COMMAND ******: " + cmd_cp_gff)
            run_cmd(cmd_cp_gff)
            run_cmd(cmd_purge_tmp)

        if self.fetch == "assembly-fna":
            print("****** NOW RUNNING COMMAND ******: " + cmd_get_assembly)
            run_cmd(cmd_get_assembly)
            print("****** NOW RUNNING COMMAND ******: " + cmd_cp_assembly)
            run_cmd(cmd_cp_assembly)
            run_cmd(cmd_purge_tmp)

        if self.fetch == "assembly-stats":
            print("****** NOW RUNNING COMMAND ******: " + cmd_get_assembly)
            run_cmd(cmd_get_assembly_stats)
            print("****** NOW RUNNING COMMAND ******: " + cmd_cp_assembly)
            run_cmd(cmd_cp_assembly_stats)
            run_cmd(cmd_purge_tmp)


        if self.fetch == "assembly-report":
            print("****** NOW RUNNING COMMAND ******: " + cmd_get_assembly)
            run_cmd(cmd_get_assembly_report)
            print("****** NOW RUNNING COMMAND ******: " + cmd_cp_assembly)
            run_cmd(cmd_cp_assembly_report)
            run_cmd(cmd_purge_tmp)

        if self.fetch == "nuclseq":
            print("****** NOW RUNNING COMMAND ******: " + cmd_get_nucleotide_fasta)
            run_cmd(cmd_get_nucleotide_fasta)

        if self.fetch == "protseq":
            print("****** NOW RUNNING COMMAND ******: " + cmd_get_nucleotide_fasta)
            run_cmd(cmd_get_protein_fasta)

class getDataNCBI(luigi.Task):
    #projectName = luigi.Parameter()
    section = luigi.Parameter(description="Choose from:['refseq','genbank']")

    fetch = luigi.Parameter(description="Choose from:['assembly-fna','assembly-gbk', 'assembly-report', 'assembly-stats', 'assembly-gff3', "
                                        "'protseq', 'nuclseq']")
    

    def requires(self):

        if self.fetch == "assembly-fna" and self.section == "refseq":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list", "ncbi_assembly_fna_refseq_accn.txt")))]]

        if self.fetch == "assembly-fna" and self.section == "genbank":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list", "ncbi_assembly_fna_genbank_accn.txt")))]]



        if self.fetch == "assembly-gff3" and self.section == "refseq":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list", "ncbi_assembly_gff3_refseq_accn.txt")))]]

        if self.fetch == "assembly-gff3" and self.section == "genbank":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list", "ncbi_assembly_gff3_genbank_accn.txt")))]]


        if self.fetch == "assembly-report" and self.section == "refseq":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list", "ncbi_assembly_report_refseq_accn.txt")))]]

        if self.fetch == "assembly-report" and self.section == "genbank":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list", "ncbi_assembly_report_genbank_accn.txt")))]]



        ####
        if self.fetch == "nuclseq" and self.section == "refseq":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(),"accession_list", "ncbi_nucl_refseq_accn.txt")))]]

        if self.fetch == "nuclseq" and self.section == "genbank":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list", "ncbi_nucl_genbank_accn.txt")))]]
        if self.fetch == "protseq" and self.section == "refseq":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list","ncbi_prot_refseq_accn.txt")))]]

        if self.fetch == "protseq" and self.section == "genbank":
            return [NCBI(fetch=self.fetch,section=self.section,
                            accession=i)
                        for i in [line.strip()
                                  for line in
                                  open((os.path.join(os.getcwd(), "accession_list","ncbi_prot_genbank_accn.txt")))]]


    def output(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        return luigi.LocalTarget('data.downlaod.complete.{t}'.format(t=timestamp))

    def run(self):
        timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
        with self.output().open('w') as outfile:
            outfile.write('download finished at {t}'.format(t=timestamp))
