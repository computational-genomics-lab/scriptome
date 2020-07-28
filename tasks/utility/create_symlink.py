#! /usr/bin/env python3
import os
import re
import shutil
from collections import OrderedDict
import optparse
from sys import exit
import psutil
import subprocess

totalcpus = psutil.cpu_count()
threads = (totalcpus-1)
mem = psutil.virtual_memory()
maxMemory= int((mem.available/1073741824) -1)


parser = optparse.OptionParser()
parser.add_option('-d', '--dataDir',
				  help="[OPTIONAL] Path to Data Directory, Default 'raw_data'",
				  dest='dataDir',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"raw_data")))
				  )

parser.add_option('-o', '--outDir',
				  help="[Optional] Path to symbolic link Data Directory, Default 'raw_data_symlink'",
				  dest='outDir',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"raw_data_symlink")))
				  )
parser.add_option('-f', '--fileExtention',
				  help="[REQUIRED] extensions for short read FASTQ files. Allowed extensions are [fq / fq.gz / fastq / fastq.gz]",
				  dest='fileExtention'
				  )

parser.add_option('-l', '--LongReadfileExtention',
				  help="extensions for a FASTQ files. Allowed extensions are [fq / fq.gz / fastq / fastq.gz]",
				  dest='LongReadfileExtention'
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
parser.add_option('-x', '--maxMemory',
				  help="[Optional Parameter] Maximum Memory in GB. Default = (available memory in GB -1)",
				  dest='maxMemory',
				  type="int",
				  default=maxMemory
				 )

options,args = parser.parse_args()
option_dict = vars(options)



dataDir = option_dict.get('dataDir')
outDir = option_dict.get('outDir')
suffix = option_dict.get('fileExtention')
lrsuffix = option_dict.get('LongReadfileExtention')




def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)


def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output



pairedendDir=os.path.abspath(os.path.join(outDir, "pe"))
matepairDir=os.path.abspath(os.path.join(outDir, "mp"))
peInterleavedDir=os.path.abspath(os.path.join(outDir, "ilpe"))
mpInterleavedDir=os.path.abspath(os.path.join(outDir, "ilmp"))
singleendDir=os.path.abspath(os.path.join(outDir, "se"))
longreadDir=os.path.abspath(os.path.join(outDir, "lr"))
email=option_dict.get('emailAddress')
port = int(option_dict.get('schedulerPort'))
cpu = int(option_dict.get('threads'))
memory = int(option_dict.get('maxMemory'))

createFolder(outDir)
createFolder(pairedendDir)
createFolder(matepairDir)
createFolder(peInterleavedDir)
createFolder(mpInterleavedDir)
createFolder(singleendDir)
createFolder(longreadDir)
createFolder("sample_list")

if os.path.isdir(dataDir):
	files = [f for f in os.listdir(dataDir) if os.path.isfile(os.path.join(dataDir, f))]
	keys = []
	fileList = re.compile(r'^(.+?).(gff|gtf|fasta|fna|fa|fasta\.gz|fna\.gz|fa\.gz|fastq|fq|fastq\.gz|fq\.gz)?$')
	for file in files:
		if fileList.search(file):
			keys.append(file)

dicts = OrderedDict ()
#keys = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

for i in keys:
        dicts[i] = str(input("Enter Data Type of {data}: \tchoose from [pe, mp, ilpe, ilmp, se, lr, gff, gtf]: ".format(data=i)))
print(dicts)

for key, val in dicts.items():
	if not os.path.exists(os.path.join(outDir, val)):
		os.mkdir(os.path.join(outDir, val))

##ln -nsf method
for key, val in dicts.items():
	dest = (os.path.join(outDir,val,key))
	src = (os.path.join(dataDir,key))
	source = os.path.abspath(src)
	destination = os.path.abspath(dest)
	escape="\'"
	print("Source:\t {skip}{src}{skip}".format(skip=escape,src=source))
	print("Destination:\t {skip}{dest}{skip}".format(skip=escape,dest=destination))
	#print("Desitnation:", '\'destination\'')

	link_cmd = "ln -nsf "
	create_symlink = "{link_cmd} {source} {destination}".format(link_cmd=link_cmd,source=source,destination=destination)
	print("****** NOW RUNNING COMMAND ******: " + create_symlink)
	print (run_cmd(create_symlink))

###########################################


if pairedendDir is not None:
	#files = [f for f in os.listdir(pairedendDir) if os.path.isfile(os.path.join(pairedendDir, f))]
	files = [f for f in os.listdir(pairedendDir) if os.path.islink(os.path.join(pairedendDir, f))]

	pe_final_list = []
	#pairedEndList = re.compile(r'^(.+?)_R[12]\.(fastq|fq|fastq\.gz|fq\.gz)?$')
	pairedEndList = re.compile(r'^(.+?)_R[12]\.(fastq\.gz|fq\.gz)?$')

	for file in files:
		if pairedEndList.search(file):
			if suffix in file:
				bname = file.split('_R')
				sample = bname[0]
				pe_final_list.append(sample)
			else:
				print("\nFile Extension Error\n")
				print("The file extention \"{suffix}\" you provided, is not matching with files at {pairedendDir}\n".format(suffix=suffix,pairedendDir=pairedendDir))
				print("List of files at {pairedendDir}".format(pairedendDir=pairedendDir))
				for file in files:
					print (file)
					exit ()

	sample_set = set(pe_final_list)
	sample_list = sorted(sample_set)
	number=len(sample_list)

	print ("\nNumber of paired-end sample(s) found:\t{number}\n".format(number=number))

	with open ((os.path.join(os.getcwd(),"sample_list",'pe_samples.lst')),'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()

#if os.path.isdir(matepairDir):
if matepairDir is not None:
	#files = [f for f in os.listdir(matepairDir) if os.path.isfile(os.path.join(matepairDir, f))]
	files = [f for f in os.listdir(matepairDir) if os.path.islink(os.path.join(matepairDir, f))]
	mp_final_list = []
	MatePairList = re.compile(r'^(.+?)_R[12]\.(fastq\.gz|fq\.gz)?$')
	#pairedEndList = re.compile(r'^(.+?)_R[12]\.(fastq|fq(?:\.gz)?)$')
	for file in files:
		if MatePairList.search(file):
			if suffix in file:
				bname = file.split('_R')
				sample = bname[0]
				mp_final_list.append(sample)
			else:
				print("\nFile Extension Error\n")
				print("The file extention \"{suffix}\" you provided, is not matching with files at {matepairDir}\n"
					  .format(suffix=suffix,matepairDir=matepairDir))
				print("List of files at {matepairDir}".format(matepairDir=matepairDir))
				for file in files:
					print (file)
					exit ()
	sample_set = set(mp_final_list)
	sample_list = sorted(sample_set)
	number=len(sample_list)

	print ("Number of mate-pair sample(s) found:\t{number}\n".format(number=number))
	with open ((os.path.join(os.getcwd(), "sample_list",'mp_samples.lst')),'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()
#if os.path.isdir(peInterleavedDir):
if peInterleavedDir is not None:
	#files = [f for f in os.listdir(peInterleavedDir) if os.path.isfile(os.path.join(peInterleavedDir, f))]
	files = [f for f in os.listdir(peInterleavedDir) if os.path.islink(os.path.join(peInterleavedDir, f))]
	pe_il_final_list = []
	#pe_il_read_list = re.compile(r'^(.+?)\.(fastq|fq|fastq\.gz|fq\.gz)?$')
	pe_il_read_list = re.compile(r'^(.+?)\.(fastq\.gz|fq\.gz)?$')
	for file in files:
		if pe_il_read_list.search(file):
			if suffix in file:
				bname = file.split('.')
				sample = bname[0]
				pe_il_final_list.append(sample)
			else:
				print("\nFile Extension Error\n")
				print("The file extention \"{suffix}\" you provided, is not matching with files at {peInterleavedDir}\n".
					format(suffix=suffix, peInterleavedDir=peInterleavedDir))
				print("List of files at {peInterleavedDir}".format(peInterleavedDir=peInterleavedDir))
				for file in files:
					print(file)
					exit()

	sample_set = set(pe_il_final_list)
	sample_list = sorted(sample_set)
	number = len(sample_list)

	print("Number of Interleaved Paired-end sample(s) found:\t{number}\n".format(number=number))

	with open((os.path.join(os.getcwd(),"sample_list", 'pe_il_samples.lst')), 'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()

#if os.path.isdir(mpInterleavedDir):
if mpInterleavedDir is not None:
	#files = [f for f in os.listdir(mpInterleavedDir) if os.path.isfile(os.path.join(mpInterleavedDir, f))]
	files = [f for f in os.listdir(mpInterleavedDir) if os.path.islink(os.path.join(mpInterleavedDir, f))]
	mp_il_final_list = []
	#mp_il_read_list = re.compile(r'^(.+?)\.(fastq|fq|fastq\.gz|fq\.gz)?$')
	mp_il_read_list = re.compile(r'^(.+?)\.(fastq\.gz|fq\.gz)?$')
	for file in files:
		if mp_il_read_list.search(file):
			if suffix in file:
				bname = file.split('.')
				sample = bname[0]
				mp_il_final_list.append(sample)
			else:
				print("\nFile Extension Error\n")
				print("The file extention \"{suffix}\" you provided, is not matching with files at {mpInterleavedDir}\n".
					format(suffix=suffix, mpInterleavedDir=mpInterleavedDir))
				print("List of files at {mpInterleavedDir}".format(mpInterleavedDir=mpInterleavedDir))
				for file in files:
					print(file)
					exit()

	sample_set = set(mp_il_final_list)
	sample_list = sorted(sample_set)
	number = len(sample_list)

	print("Number of Interleaved Mate-pair sample(s) found:\t{number}\n".format(number=number))

	with open((os.path.join(os.getcwd(),"sample_list", 'mp_il_samples.lst')), 'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()

##########
if singleendDir is not None:
	#files = [f for f in os.listdir(singleendDir) if os.path.isfile(os.path.join(singleendDir, f))]
	files = [f for f in os.listdir(singleendDir) if os.path.islink(os.path.join(singleendDir, f))]
	se_final_list = []

	#se_read_list = re.compile(r'^(.+?)\.(fastq|fq|fastq\.gz|fq\.gz)?$')
	se_read_list = re.compile(r'^(.+?)\.(fastq\.gz|fq\.gz)?$')
	for file in files:
		if mp_il_read_list.search(file):
			if suffix in file:
				bname = file.split('.')
				sample = bname[0]
				se_final_list.append(sample)
			else:
				print("\nFile Extension Error\n")
				print("The file extention \"{suffix}\" you provided is not matching with files at {singleendDir}\n".
					format(suffix=suffix, singleendDir=singleendDir))
				print("List of files at {singleendDir}".format(singleendDir=singleendDir))
				for file in files:
					print(file)
					exit()

	sample_set = set(se_final_list)
	sample_list = sorted(sample_set)
	number = len(sample_list)

	print("Number of Single-end sample(s) found:\t{number}\n".format(number=number))

	with open((os.path.join(os.getcwd(),"sample_list",'se_samples.lst')), 'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()
##########
if longreadDir is not None:
	#files = [f for f in os.listdir(longreadDir) if os.path.isfile(os.path.join(longreadDir, f))]
	files = [f for f in os.listdir(longreadDir) if os.path.islink(os.path.join(longreadDir, f))]
	lr_final_list = []

	lr_read_list = re.compile(r'^(.+?)\.(fasta\.gz|fastq\.gz|fq\.gz)?$')
	for file in files:
		if lr_read_list.search(file):
			if lrsuffix in file:
				bname = file.split('.')
				sample = bname[0]
				lr_final_list.append(sample)
			else:
				print("\nFile Extension Error\n")
				print("The file extention \"{lrsuffix}\" you provided, is not matching with files at {longreadDir}\n".
					format(lrsuffix=lrsuffix, longreadDir=longreadDir))
				print("List of files at {longreadDir}".format(longreadDir=longreadDir))
				for file in files:
					print(file)
					exit()

	sample_set = set(lr_final_list)
	sample_list = sorted(sample_set)
	number = len(sample_list)

	print("Number of long read sample(s) found:\t{number}\n".format(number=number))

	with open((os.path.join(os.getcwd(),"sample_list", 'lr_samples.lst')), 'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()
##########


with open('luigi.cfg', 'w') as fh:
	fh.write('[core]\n')
	fh.write('default-scheduler-port:{port}\n'.format(port=port))
	fh.write('error-email={email}\n\n'.format(email=email))
	fh.write('[GlobalParameter]\n')
	fh.write('singleendDir={singleendDir}\n'.format(singleendDir=singleendDir))
	fh.write('pairedendDir={pairedendDir}\n'.format(pairedendDir=pairedendDir))
	fh.write('matepairDir={matepairDir}\n'.format(matepairDir=matepairDir))
	fh.write('mpInterleavedDir={mpInterleavedDir}\n'.format(mpInterleavedDir=mpInterleavedDir))
	fh.write('peInterleavedDir={peInterleavedDir}\n'.format(peInterleavedDir=peInterleavedDir))
	fh.write('suffix={suffix}\n'.format(suffix=suffix))
	fh.write('longreadDir={longreadDir}\n'.format(longreadDir=longreadDir))
	fh.write('lrsuffix={lrsuffix}\n'.format(lrsuffix=lrsuffix))
	#fh.write('genomeDir={genomeDir}\n'.format(genomeDir=genomeDir))
	#fh.write('projectName={projectName}\n'.format(projectName=projectName))
	fh.write('threads={cpu}\n'.format(cpu=cpu))
	fh.write('maxMemory={memory}\n'.format(memory=memory))
	fh.close()
