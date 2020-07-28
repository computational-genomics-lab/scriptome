import os
import luigi
import os
import time
import subprocess
import pandas as pd
from luigi import Parameter

from tasks.assembly.kmergenie import kmergenie_formater_bbduk
from tasks.assembly.kmergenie import kmergenie_formater_reformat
from tasks.assembly.kmergenie import optimal_kmer
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import filtlong
from tasks.readCleaning.reFormatReads import reformat
suffix='.fastq'
old_left_read_name_suffix = '_R1.fastq'
old_rigt_read_name_suffix = '_R2.fastq'
new_left_read_name_suffix = '_1.fastq'
new_rigt_read_name_suffix = '_2.fastq'

cleaned_read_folder = os.path.join(os.getcwd(), "CleanedReads", "Cleaned_PE_Reads" + "/")
verified_read_folder = os.path.join(os.getcwd(), "VerifiedReads", "Verified_PE_Reads" + "/")
symlink_read_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "input_reads" + "/")
genome_bin_folder =  os.path.join(os.getcwd(), "MetagenomeAnalysis", "binning" + "/")

task_log_folder =  os.path.join(os.getcwd(), "MetagenomeAnalysis", "task_logs" + "/")


def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)
createFolder(symlink_read_folder)
createFolder(task_log_folder)

def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
						 shell=True,
						 universal_newlines=True,
						 stdout=subprocess.PIPE,
						 executable='/bin/bash')
	output = p.communicate()[0]
	return output
def clean_read_symlink(input_file):
	with open(input_file) as ifh:
		sample_name_list = ifh.read().splitlines()
		dir = os.listdir(symlink_read_folder)
		if len(dir) != 0: 
				print("symolic link directory not empty" )
				cmd_rm_symlink_dir="rm -rf {symlink_read_folder}".format(symlink_read_folder=symlink_read_folder) 
				print("****** NOW RUNNING COMMAND ******: " + cmd_rm_symlink_dir)
				print (run_cmd(cmd_rm_symlink_dir))
				createFolder(symlink_read_folder)		

		for sample in sample_name_list:	
			if not (len(dir) ==  len(sample_name_list)): 
				os.symlink(os.path.join(cleaned_read_folder,sample + old_left_read_name_suffix),os.path.join(symlink_read_folder,sample+new_left_read_name_suffix))
				os.symlink(os.path.join(cleaned_read_folder,sample + old_rigt_read_name_suffix),os.path.join(symlink_read_folder,sample+new_rigt_read_name_suffix))	
		print('symlinks created successfully')
		cmd_tree="tree {symlink_read_folder}".format(symlink_read_folder=symlink_read_folder)
		print (run_cmd(cmd_tree))

def verified_read_symlink(input_file):
	with open(input_file) as ifh:
		sample_name_list = ifh.read().splitlines()
		dir = os.listdir(symlink_read_folder)
		if len(dir) != 0: 
				print("symolic link directory not empty" )
				cmd_rm_symlink_dir="rm -rf {symlink_read_folder}".format(symlink_read_folder=symlink_read_folder) 
				print("****** NOW RUNNING COMMAND ******: " + cmd_rm_symlink_dir)
				print (run_cmd(cmd_rm_symlink_dir))
				createFolder(symlink_read_folder)		

		for sample in sample_name_list:	
			if not (len(dir) ==  len(sample_name_list)): 
				os.symlink(os.path.join(verified_read_folder,sample + old_left_read_name_suffix),os.path.join(symlink_read_folder,sample+new_left_read_name_suffix))
				os.symlink(os.path.join(verified_read_folder,sample + old_rigt_read_name_suffix),os.path.join(symlink_read_folder,sample+new_rigt_read_name_suffix))	
		print('symlinks created successfully')
		cmd_tree="tree {symlink_read_folder}".format(symlink_read_folder=symlink_read_folder)
		print (run_cmd(cmd_tree))



############################################
class GlobalParameter(luigi.Config):
	read_library_type=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	adapter=luigi.Parameter()

class metaphlan(luigi.Task):
	project_name=luigi.Parameter(default="MetagenomeAnalysis")
	#assembly_type=luigi.ChoiceParameter(default="single",choices=["single", "co"], var_type=str)
	adapter = GlobalParameter().adapter
	threads = GlobalParameter().threads
	max_memory = GlobalParameter().maxMemory
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	read_library_type = GlobalParameter().read_library_type


	def requires(self):
		if self.pre_process_reads=="yes":
			return [bbduk(read_library_type="pe",
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		
		if self.pre_process_reads=="no":
			return [reformat(read_library_type="pe",
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
		

	def output(self):
		inDir = os.path.join(os.getcwd(), "MetagenomeAnalysis", "input_reads" + "/")
		readlist=os.listdir(inDir)

		for read in readlist:	
			if read.endswith("_1.fastq"):
				filename,file_extn=read.split(".",1)[0],read.split('.',1)[1]
				forward_read=read
				reverse_read=read.replace("_1.fastq", "_2.fastq")
				sample=read.split("_1",1)[0]

				metaphlan_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" )
				return {'out1': luigi.LocalTarget(metaphlan_folder + "/" + "profiled_samples" + "/" + sample+".txt"),
						'out2': luigi.LocalTarget(metaphlan_folder + "/" + "profiled_samples" + "/" + sample+".out")}
		
	def run(self):
		metaphlan_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" )
		input_sample_list=os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")


		if self.pre_process_reads == "yes":
			#cmd_create_clean_read_symlink = clean_read_symlink(input_sample_list)
			clean_read_symlink(input_sample_list)

		if self.pre_process_reads == "no":
			#cmd_create_verified_read_symlink = verified_read_symlink(input_sample_list)
			verified_read_symlink(input_sample_list)

		inDir=os.path.join(os.getcwd(), "MetagenomeAnalysis", "input_reads" + "/")
		readlist=os.listdir(inDir)


		for read in readlist:
			allowed_extn=["_1.fastq"]
			if any (read.endswith(ext) for ext in allowed_extn):
				filename,file_extn=read.split(".",1)[0],read.split('.',1)[1]
				forward_read=read
				reverse_read=read.replace("_1.fastq", "_2.fastq")

				sample_name=read.split("_1",1)[0]
				print(sample_name)
				print(forward_read)
				print(reverse_read)
				print(file_extn)

				run_metaphlan="[ -d {metaphlan_folder}/profiled_samples ] || mkdir -p {metaphlan_folder}/profiled_samples ; " \
					  "metaphlan --input_type fastq " \
					  "{inDir}{forward_read},{inDir}{reverse_read} " \
					  "--nproc {threads} " \
					  "--bt2_ps very-sensitive " \
					  "--bowtie2out {metaphlan_folder}/profiled_samples/{sample}.out 2>&1 | tee {metaphlan_folder}/profiled_samples/{sample}.txt ".format(
						inDir=inDir,metaphlan_folder=metaphlan_folder,sample=sample_name,reverse_read=reverse_read,forward_read=forward_read,threads=self.threads)

				print("****NOW RUNNING COMMAND****:" + run_metaphlan)
				print(run_cmd(run_metaphlan))

class graphlan(luigi.Task):
	project_name=luigi.Parameter(default="MetagenomeAnalysis")
	#assembly_type=luigi.ChoiceParameter(default="single",choices=["single", "co"], var_type=str)
	adapter = GlobalParameter().adapter
	threads = GlobalParameter().threads
	max_memory = GlobalParameter().maxMemory
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	read_library_type = GlobalParameter().read_library_type


	def requires(self):
		return [metaphlan(pre_process_reads=self.pre_process_reads)]

	def output(self):
		metaphlan_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"profiled_samples" + "/")
		outDir=os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"data_for_images" + "/")
		graphlan_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"figures" + "/")

		return {'out1': luigi.LocalTarget(outDir + "/" + "merged_abundance_table.txt"),
				'out2': luigi.LocalTarget(outDir + "/" + "otu_table_ampvis2.txt"),
				'out3': luigi.LocalTarget(outDir + "/" + "otu_table_phyloseq.txt"),
				'out4': luigi.LocalTarget(outDir + "/" + "otu_table_phyloseq.biom"),
				'out5': luigi.LocalTarget(graphlan_folder + "/" + "hclust_abundance_heatmap_species.png"),
				'out6': luigi.LocalTarget(graphlan_folder + "/" + "graphlan_merged_abundance.pdf")}

	def run(self):
		metaphlan_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"profiled_samples" + "/")
		outDir=os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"data_for_images" + "/")
		graphlan_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"figures" + "/")


		run_merge_metaphln_tables="[ -d {outDir} ] || mkdir -p {outDir} ; " \
								  "merge_metaphlan_tables.py {metaphlan_folder}*.txt > {outDir}merged_abundance_table.txt".format(outDir=outDir,metaphlan_folder=metaphlan_folder)
		
		print("****NOW RUNNING COMMAND****:" + run_merge_metaphln_tables)
		print(run_cmd(run_merge_metaphln_tables))


		cmd_prepare_file="grep -E '(s__)' {outDir}/merged_abundance_table.txt |cut -f2 --complement | grep -v 't__' | sed 's/^.*s__//g' " \
						 "> {outDir}/merged_abundance_table_species.txt".format(outDir=outDir)

		cmd_generate_header="cat {outDir}/merged_abundance_table.txt | head -n 2 | tail -n 1 | cut -f1 --complement " \
							"> {outDir}/header.txt".format(outDir=outDir)

		cmd_merge_header="cat {outDir}/header.txt {outDir}/merged_abundance_table_species.txt " \
						 "> {outDir}/abundance_table_species.txt".format(outDir=outDir)

		cmd_taxonomy_phyloseq="cat {outDir}/merged_abundance_table.txt | grep -E 's__|clade' " \
							"| grep -v '#'| grep -v 'MetaPhlAn'| grep -v 'WARNING'| " \
							"sed -e 's/clade_name/taxonomy/'|sed -e 's/|/;/g'|cut -f1 " \
							 "> {outDir}/taxonomy_phyloseq.txt".format(outDir=outDir)

		cmd_sample_abundance="cat {outDir}/merged_abundance_table.txt | grep -E 's__|clade' " \
							 "| grep -v '#'| grep -v 'MetaPhlAn'| grep -v 'WARNING'| " \
							 "sed -e 's/clade_name/taxonomy/'|sed -e 's/|/;/g'|cut -f1,2 --complement " \
							 "> {outDir}/sample_abundance.txt".format(outDir=outDir)

		cmd_otu_id_phyloseq="cat {outDir}/taxonomy_phyloseq.txt |sed 's/^.*s__//g' | sed -e 's/taxonomy/#OTU ID/' " \
							"> {outDir}/otu_id_phyloseq.txt".format(outDir=outDir)

		cmd_otu_table_pgyloseq="paste {outDir}/otu_id_phyloseq.txt {outDir}/sample_abundance.txt {outDir}/taxonomy_phyloseq.txt " \
							   "> {outDir}otu_table_phyloseq.txt".format(outDir=outDir)

		cmd_phyloseq_biom="biom convert -i  {outDir}/otu_table_phyloseq.txt " \
						  "-o  {outDir}/otu_table_phyloseq.biom --table-type='OTU table' " \
						  "--process-obs-metadata taxonomy --to-json".format(outDir=outDir)

		cmd_taxonomy_ampvis2="cat {outDir}/taxonomy_phyloseq.txt |sed -e 's/;/;\t/g' | " \
							 "sed -e '1s/taxonomy/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/g' " \
							"> {outDir}/taxonomy_ampvis2.txt".format(outDir=outDir)

		cmd_otu_phyloseq="cat {outDir}/taxonomy_phyloseq.txt |sed 's/^.*s__//g' " \
						 "| sed -e 's/taxonomy/OTU/' > {outDir}/otu_id_ampvis2.txt".format(outDir=outDir)

		cmd_otu_ampvis2="paste {outDir}/otu_id_ampvis2.txt {outDir}/sample_abundance.txt {outDir}/taxonomy_ampvis2.txt" \
						" > {outDir}/otu_table_ampvis2.txt".format(outDir=outDir)

		
		print("****NOW RUNNING COMMAND****:" + cmd_taxonomy_phyloseq)
		print(run_cmd(cmd_taxonomy_phyloseq))		
		print("****NOW RUNNING COMMAND****:" + cmd_sample_abundance)
		print(run_cmd(cmd_sample_abundance))
		print("****NOW RUNNING COMMAND****:" + cmd_otu_id_phyloseq)
		print(run_cmd(cmd_otu_id_phyloseq))	
		print("****NOW RUNNING COMMAND****:" + cmd_otu_table_pgyloseq)
		print(run_cmd(cmd_otu_table_pgyloseq))	
		print("****NOW RUNNING COMMAND****:" + cmd_phyloseq_biom)
		print(run_cmd(cmd_phyloseq_biom))	
		print("****NOW RUNNING COMMAND****:" + cmd_taxonomy_ampvis2)
		print(run_cmd(cmd_taxonomy_ampvis2))	
		print("****NOW RUNNING COMMAND****:" + cmd_otu_phyloseq)
		print(run_cmd(cmd_otu_phyloseq))	
		print("****NOW RUNNING COMMAND****:" + cmd_otu_ampvis2)
		print(run_cmd(cmd_otu_ampvis2))	

		print("****NOW RUNNING COMMAND****:" + cmd_prepare_file)
		print(run_cmd(cmd_prepare_file))

		print("****NOW RUNNING COMMAND****:" + cmd_generate_header)
		print(run_cmd(cmd_generate_header))


		print("****NOW RUNNING COMMAND****:" + cmd_merge_header)
		print(run_cmd(cmd_merge_header))



		cmd_run_hclust2="[ -d {graphlan_folder} ] || mkdir -p {graphlan_folder} ; " \
					"hclust2.py -i {outDir}/abundance_table_species.txt " \
					"-o {graphlan_folder}/hclust_abundance_heatmap_species.png " \
					"--ftop 25 " \
					"--f_dist_f braycurtis " \
					"--s_dist_f braycurtis " \
					"--cell_aspect_ratio 0.5 -l " \
					"--flabel_size 6 " \
					"--slabel_size 6 " \
					"--max_flabel_len 100 " \
					"--max_slabel_len 100 --minv 0.1 --dpi 300 ".format(outDir=outDir,graphlan_folder=graphlan_folder)

	

		cmd_abundance_reformat="tail -n +5 {outDir}/merged_abundance_table.txt | cut -f1,3- " \
							   "> {outDir}/merged_abundance_table_reformatted.txt".format(outDir=outDir)
		
		print("****NOW RUNNING COMMAND****:" + cmd_abundance_reformat)
		print(run_cmd(cmd_abundance_reformat))

		cmd_run_graphlan_prep="export2graphlan.py --skip_rows 1 " \
					"-i {outDir}/merged_abundance_table_reformatted.txt " \
					"--tree {outDir}/merged_abundance.tree.txt " \
					"--annotation {outDir}/merged_abundance.annot.txt " \
					"--most_abundant 100 " \
					"--abundance_threshold 1 " \
					"--least_biomarkers 10 " \
					"--annotations 5,6 " \
					"--external_annotations 7 "  \
					"--min_clade_size 1 ".format(outDir=outDir)


		cmd_run_graphlan_annotate="[ -d {outDir} ] || mkdir -p {outDir} ; " \
					"graphlan_annotate.py " \
					"--annot {outDir}/merged_abundance.annot.txt " \
					"{outDir}/merged_abundance.tree.txt " \
					"{outDir}/merged_abundance.xml ".format(outDir=outDir)


		cmd_run_graphlan="[ -d {graphlan_folder} ] || mkdir -p {graphlan_folder} ; " \
					"graphlan.py " \
					"--dpi 300 " \
					"--size 15 " \
					"--format pdf " \
					"{outDir}/merged_abundance.xml " \
					"{graphlan_folder}/graphlan_merged_abundance.pdf ".format(outDir=outDir,graphlan_folder=graphlan_folder)

		
		
		print("****NOW RUNNING COMMAND****:" + cmd_run_hclust2)
		print(run_cmd(cmd_run_hclust2))

		print("****NOW RUNNING COMMAND****:" + cmd_run_graphlan_prep)
		print(run_cmd(cmd_run_graphlan_prep))

		print("****NOW RUNNING COMMAND****:" + cmd_run_graphlan_annotate)
		print(run_cmd(cmd_run_graphlan_annotate))

		print("****NOW RUNNING COMMAND****:" + cmd_run_graphlan)
		print(run_cmd(cmd_run_graphlan))


class profileTaxonomy(luigi.Task):
	project_name=luigi.Parameter(default="MetagenomeAnalysis")
	adapter = GlobalParameter().adapter
	threads = GlobalParameter().threads
	max_memory = GlobalParameter().maxMemory
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	read_library_type = GlobalParameter().read_library_type
	condition_column=luigi.Parameter()
	alpha=luigi.FloatParameter(default=0.05)


	def requires(self):
		return [graphlan(pre_process_reads=self.pre_process_reads)]

	def output(self):
		phyloseq_image_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"figures" + "/")
		return {'out1': luigi.LocalTarget(phyloseq_image_folder + "/" + "volcano_plot_FC_2_P_05.tiff.tiff")}


	def run(self):
		phyloseq_folder = os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" + "/")
		inDir=os.path.join(os.getcwd(), "MetagenomeAnalysis", "metaphlan_analysis" ,"data_for_images" + "/")
		map_file=os.path.join(os.getcwd(), "sample_list","map.txt")

		cmd_run_phyloseq="[ -d {phyloseq_folder} ] || mkdir -p {phyloseq_folder} ;  cd {phyloseq_folder} ;" \
						"phyloseq.r -t {map_file} " \
						"-P {inDir}/otu_table_phyloseq.biom " \
						"-v {condition_column} " \
						"-a {alpha} " \
						"-A {inDir}/otu_table_ampvis2.txt".format(map_file=map_file,
																 inDir=inDir,
																 alpha=self.alpha,
																 phyloseq_folder=phyloseq_folder,
																 condition_column=self.condition_column)
		print("****NOW RUNNING COMMAND****:" + cmd_run_phyloseq)
		print(run_cmd(cmd_run_phyloseq))