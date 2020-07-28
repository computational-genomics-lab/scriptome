#!/usr/bin/env python
import optparse
import re
import os
from sys import exit
import pandas as pd
import psutil
from collections import OrderedDict

########################################################################################################################

required="fileExtention domain  projectName" \
		 "".split()
parser = optparse.OptionParser()

totalcpus = psutil.cpu_count()
threads = (totalcpus-1)

mem = psutil.virtual_memory()
maxMemory= int((mem.available/1073741824) -1)


parser.add_option('-P', '--pairedendDir',
				  help="Path to Directory containing short insert paired-end reads",
				  dest='pairedendDir'
				  )
parser.add_option('-M', '--matepairDir',
				  help="Path to Directory containing long insert matepair reads",
				  dest='matepairDir'
				  )
parser.add_option('-L', '--longreadDir',
				  help="Path to Directory containing long nanopore reads",
				  dest='longreadDir'
				  )
parser.add_option('-T', '--trustedContigs',
				  help="Path to Directory containing trusted contigs from any other assemblers",
				  dest='trustedContigs'
				  )
parser.add_option('-G', '--genomeDir',
				  help="Path to Directory containing reference genome (.fna)",
				  dest='genomeDir'
				 )
parser.add_option('-x', '--fileExtention',
				  help="extensions for a FASTQ files. Allowed extensions are [fq / fq.gz / fastq / fastq.gz]",
				  dest='fileExtention'
				  )
parser.add_option('-d', '--domain',
				  help="Organism Domain =string[prokaryote OR eukaryote]",
				  dest='domain'
				  )
parser.add_option('-g', '--approxGenomeSize',
				  help="Approximate Genome Size in MB =float[10.7]",
				  dest='approxGenomeSize'
				  )
parser.add_option('-z', '--highlyHeterozygous',
				  help="whether highly heterozygous genome?  =str[yes / no]",
				  dest='highlyHeterozygous',
				  default='no'
				  )
parser.add_option('-a', '--adapterFile',
				  help="Path to adapter file",
				  dest='adapterFile'
				  )
parser.add_option('-o', '--projectName',
				  help="Path to the Project Directory ",
				  dest='projectName',
				  default='MyProject'
				 )
parser.add_option('-s', '--schedulerPort',
				  help="[Optional Parameter] Scheduler Port Number. default =int[8888] ",
				  type="int",
				  default=8082
				 )
parser.add_option('-e', '--emailAddress',
				  help="Provide your email address =email[abc@xyz.com]",
				  dest='emailAddress'
				 )
parser.add_option('-t', '--threads',
				  help="[Optional Parameter, ] Number of threads. Default = (total threads -1)",
				  dest='threads',
				  type="int",
				  default = threads)
parser.add_option('-m', '--maxMemory',
				  help="[Optional Parameter] Maximum Memory in GB. Default = (available memory in GB -1)",
				  dest='maxMemory',
				  type="int",
				  default=maxMemory
				 )
options,args = parser.parse_args()
for r in required:
	if options.__dict__[r] is None:
		parser.error("parameter %s required" %r)


option_dict = vars(options)
pairedendDir = option_dict.get('pairedendDir')
matepairDir = option_dict.get('matepairDir')
longreadDir = option_dict.get('longreadDir')
adapterPath =option_dict.get('adapterFile')
trustedContigsDir = option_dict.get('trustedContigs')
suffix = option_dict.get('fileExtention')
domain = option_dict.get('domain')
genomeSize = float(option_dict.get('approxGenomeSize'))
heterozygocity = option_dict.get('highlyHeterozygous')
genomeDir = option_dict.get('genomeDir')
projectName = option_dict.get('projectName')
email=option_dict.get('emailAddress')
port = int(option_dict.get('schedulerPort'))
cpu = int(option_dict.get('threads'))
memory = int(option_dict.get('maxMemory'))

if domain == "prokaryote":
	if None not in (pairedendDir, matepairDir, trustedContigsDir, longreadDir):
		assembler = "spades"

	if all(value is None for value in [pairedendDir, matepairDir, trustedContigsDir]):
		assembler = "flye"

	if (all(value is None for value in [matepairDir, trustedContigsDir,longreadDir])):
		assembler = "skesa" 

	if ((None not in (pairedendDir, longreadDir)) and None in (matepairDir, trustedContigsDir)):
		assembler = "unicycler"
	
	if ((None not in (pairedendDir, matepairDir)) and (heterozygocity == "yes")):
		assembler = "platanusb"


if domain == "eukaryote":
	if (((None not in (pairedendDir, matepairDir, longreadDir)) or None not in (pairedendDir, matepairDir)) and (genomeSize <= 80)):
		assembler = "spades"

	if (((None not in (pairedendDir, matepairDir, longreadDir)) or None not in (pairedendDir, matepairDir)) and (genomeSize > 80)):
		assembler = "abyss"
	
	if (((None not in (pairedendDir, longreadDir)) and None in (trustedContigsDir, matepairDir)) and (genomeSize > 80)):
		assembler = "dbg2olc"

	if (((None not in (pairedendDir, matepairDir, longreadDir)) or None not in (pairedendDir, matepairDir)) and ((genomeSize > 80) and (heterozygocity == "yes"))):
		assembler = "platanus"
	
	if all(value is None for value in [pairedendDir, matepairDir, trustedContigsDir]):
		assembler = "flye"

	
def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)

# Example
createFolder(projectName)


with open('luigi.cfg', 'w') as fh:
	fh.write('[core]\n')
	fh.write('default-scheduler-port:{port}\n'.format(port=port))
	fh.write('error-email={email}\n\n'.format(email=email))
	fh.write('[GlobalParameter]\n')
	fh.write('pairedendDir={pairedendDir}\n'.format(pairedendDir=pairedendDir))
	fh.write('matepairDir={matepairDir}\n'.format(matepairDir=matepairDir))
	fh.write('longreadDir={longreadDir}\n'.format(longreadDir=longreadDir))
	fh.write('trustedContigsDir={trustedContigsDir}\n'.format(trustedContigsDir=trustedContigsDir))
	fh.write('suffix={suffix}\n'.format(suffix=suffix))
	fh.write('domain={domain}\n'.format(domain=domain))
	fh.write('assembler={assembler}\n'.format(assembler=assembler))
	fh.write('heterozygous={heterozygocity}\n'.format(heterozygocity=heterozygocity))
	fh.write('adapter={adapterPath}\n'.format(adapterPath=adapterPath))
	fh.write('genomeDir={genomeDir}\n'.format(genomeDir=genomeDir))
	fh.write('projectName={projectName}\n'.format(projectName=projectName))
	fh.write('threads={cpu}\n'.format(cpu=cpu))
	fh.write('maxMemory={memory}\n'.format(memory=memory))
	fh.close()

if not [ pairedendDir == None ]:
	if os.path.isdir(pairedendDir):
		files = [f for f in os.listdir(pairedendDir) if os.path.isfile(os.path.join(pairedendDir, f))]
		pe_final_list = set()
		peList = re.compile(r'^(.+?)_R[12]\.(fastq|fq|fastq\.gz|fq\.gz)?$')
		for file in files:
			if peList.search(file):
				if suffix in file:
					bname = file.split('_R')
					sample = bname[0]
					pe_final_list.add(sample)
			
				else:
					print("\nFile Extension Error\n")
					print("The file extention \"{suffix}\" you provided with -x argument, is not matching with files at {pairedendDir}\n".format(suffix=suffix,pairedendDir=pairedendDir))
					print("List of files at {pairedendDir}".format(pairedendDir=pairedendDir))
					for file in files:
						print (file)
					exit ()

	with open ((os.path.join(os.getcwd(),projectName,'pe_samples.lst')),'w') as file:
		for sample in pe_final_list:
			file.write("%s\n" % sample)
	file.close()

if not [ matepairDir == None ]:
	if os.path.isdir(matepairDir):
		files = [f for f in os.listdir(matepairDir) if os.path.isfile(os.path.join(matepairDir, f))]
		mp_final_list = set()
		mpList = re.compile(r'^(.+?)_R[12]\.(fastq|fq|fastq\.gz|fq\.gz)?$')
		for file in files:
			if mpList.search(file):
				if suffix in file:
					bname = file.split('_R')
					sample = bname[0]
					mp_final_list.add(sample)
				#for file in final_list:
					#print (file)
				else:
					print("\nFile Extension Error\n")
					print("The file extention \"{suffix}\" you provided with -x argument, is not matching with files at {matepairDir}\n".format(suffix=suffix,matepairDir=matepairDir))
					print("List of files at {matepairDir}".format(matepairDir=matepairDir))
					for file in files:
						print (file)
					exit ()

	with open ((os.path.join(os.getcwd(),projectName,'mp_samples.lst')),'w') as file:
		for sample in mp_final_list:
			file.write("%s\n" % sample)
	file.close()

if not [ longreadDir == None ]:
	if os.path.isdir(longreadDir):
		files = [f for f in os.listdir(longreadDir) if os.path.isfile(os.path.join(longreadDir, f))]
	
	with open ((os.path.join(os.getcwd(),projectName,'lr_sample.lst')),'w') as file:
		for sample in files:
			file.write("%s\n" % sample)
	file.close()


if not [ longreadDir == None ]:
	if os.path.isdir(longreadDir):
		files = [f for f in os.listdir(longreadDir) if os.path.isfile(os.path.join(longreadDir, f))]
	
	with open ((os.path.join(os.getcwd(),projectName,'lr_sample.lst')),'w') as file:
		for sample in files:
			file.write("%s\n" % sample)
	file.close()