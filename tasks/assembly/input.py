
import os

def spades_pe_formater(samplefile):
	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		left_read_name_prefix = '--pe1-1 '
		right_read_name_suffix = '_R2.fastq'
		right_read_name_prefix = '--pe1-2 '


		pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")

		left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
		left_read_name_list = [left_read_name_prefix + x for x in left_read_name_list]

		right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
		right_read_name_list = [right_read_name_prefix + x for x in right_read_name_list]

		result = [item for sublist in zip(left_read_name_list, right_read_name_list) for item in sublist]

		parse_string = " ".join(result)

		return parse_string

'''
infile = "./pe_samples.list"

spades_input = spades_pe_formater(infile)
print ("SPADES Assembler --Input Read Format for paired-end reads only")
print ("----------------------------------------------------------------------------------------")
print (spades_input)
print ("\n")
print ("\n")
'''

def spades_pe_mp_formater(pefile,mpfile):
	pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")
	mp_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_MP_Reads" + "/")

	with open(pefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		left_read_name_prefix = '--pe1-1 '
		right_read_name_suffix = '_R2.fastq'
		right_read_name_prefix = '--pe1-2 '

		left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
		left_read_name_list = [left_read_name_prefix + x for x in left_read_name_list]

		right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
		right_read_name_list = [right_read_name_prefix + x for x in right_read_name_list]
		pe_result = [item for sublist in zip(left_read_name_list, right_read_name_list) for item in sublist]


	with open(mpfile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		left_read_name_prefix = '--mp1-1 '
		right_read_name_suffix = '_R2.fastq'
		right_read_name_prefix = '--mp1-2 '

		left_read_name_list = [mp_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
		left_read_name_list = [left_read_name_prefix + x for x in left_read_name_list]

		right_read_name_list = [mp_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
		right_read_name_list = [right_read_name_prefix + x for x in right_read_name_list]
		mp_result = [item for sublist in zip(left_read_name_list, right_read_name_list) for item in sublist]

		

		pe_parse_string = " ".join(pe_result)
		mp_parse_string = " ".join(mp_result)

		parse_string = pe_parse_string + " " +mp_parse_string


		return parse_string

'''
pefile = "./pe_samples.list"
mpfile = "./pe_samples.list"

spades_input = spades_pe_mp_formater(pefile,mpfile)
print ("SPADES Assembler --Input Read Format for Paired-end and mate-pair reads as input ")
print ("----------------------------------------------------------------------------------------")
print (spades_input)
print ("\n")
print ("\n")
'''

def unicycler_pe_formater(samplefile):
	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		left_read_name_prefix = '--short1 '
		right_read_name_suffix = '_R2.fastq'
		right_read_name_prefix = '--short2 '


		pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")

		left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
		left_read_name_list = [left_read_name_prefix + x for x in left_read_name_list]

		right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
		right_read_name_list = [right_read_name_prefix + x for x in right_read_name_list]

		result = [item for sublist in zip(left_read_name_list, right_read_name_list) for item in sublist]

		parse_string = " ".join(result)

		return parse_string

'''
infile = "./pe_samples.list"

spades_input = unicycler_pe_formater(infile)
print ("Unicycler Assembler --Input Read Format for short paired-end reads only")
print ("----------------------------------------------------------------------------------------")
print (spades_input)
print ("\n")
print ("\n")
'''



def unicycler_pe_lr_formater(pefile,sefile):
	pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")
	lr_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_Long_Reads" + "/")

	with open(pefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		left_read_name_prefix = '--short1 '
		right_read_name_suffix = '_R2.fastq'
		right_read_name_prefix = '--short2 '

		left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
		left_read_name_list = [left_read_name_prefix + x for x in left_read_name_list]

		right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
		right_read_name_list = [right_read_name_prefix + x for x in right_read_name_list]
		parse_string = [item for sublist in zip(left_read_name_list, right_read_name_list) for item in sublist]
		pe_parse_string = ','.join(parse_string)


	with open(sefile) as fh:
		sample_name_list = fh.read().splitlines()
		read_name_suffix = '.fastq'
		

		read_name_list = [ x + read_name_suffix for x in sample_name_list]
		
		reads = ','.join(read_name_list)
		
		lr_parse_string = "--unpaired " + lr_cleaned_read_folder + reads

		parse_string = pe_parse_string + " " + lr_parse_string


		return parse_string


'''
pefile = "./pe_samples.list"
sefile = "./pe_samples.list"

unicycler_pe_se_input = unicycler_pe_se_formater(pefile,sefile)
print ("Input Read Format for Unicycler PE SE")
print ("----------------------------------------------------------------------------------------")
print (unicycler_pe_se_input)
print ("\n")
print ("\n")
'''



def rockhopper_formater(samplefile):
	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'

		left_read_name_list = [x + left_read_name_suffix for x in sample_name_list]
		right_read_name_list = [x + right_read_name_suffix for x in sample_name_list]

		result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]

		result1 = [x+"%" +y for x, y in result]
		parse_string = ','.join(result1)
		return parse_string

'''
infile = "./pe_samples.list"

rockhopper_input = rockhopper_formater(infile)

print ("Input Read Format for RockHopper")
print ("----------------------------------------------------------------------------------------")
print (rockhopper_input)
print ("\n")
print ("\n")

'''

def skesa_formater(samplefile):
	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'

		left_read_name_list = [x + left_read_name_suffix for x in sample_name_list]
		right_read_name_list = [x + right_read_name_suffix for x in sample_name_list]

		pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")


		result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]

		result1 = ["--fastq " + pe_cleaned_read_folder + x +"," + pe_cleaned_read_folder +y for x, y in result]
		parse_string = ' '.join(result1)
		return parse_string


'''
infile = "./pe_samples.list"

skesa_input = skesa_formater(infile)

print ("Input Read Format for SKESA Assembler")
print ("----------------------------------------------------------------------------------------")
print (skesa_input)
print ("\n")
print ("\n")
'''

def abyss_pe(samplefile):
	with open(samplefile) as fh:
		pe_sample_name_list = fh.read().splitlines()
		pe_samples=[x for x in pe_sample_name_list]

		pe_lib_string ="lib=" + "'" + ' '.join(pe_samples) + "'"

		pe_left_read_name_suffix = '_R1.fastq'
		pe_right_read_name_suffix = '_R2.fastq'

		pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
		pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
		pe_samples=[x for x in pe_sample_name_list]

		pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]

		pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")


		pe_result1 =   [' '+ x +'=' + "'"+ pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for x, y,z in pe_result]

		#lib_name = ' '.join(samples) + "="

		pe_result2 = ' '.join(pe_result1)
		
		pe_parse_string = pe_lib_string  + pe_result2
		return pe_parse_string

'''
infile = "./pe_samples.list"

abyss_input = abyss_pe(infile)
print ("Input PE Libraries for ABYSS2 Assembler")
print ("----------------------------------------------------------------------------------------")
print (abyss_input)
print ("\n")
print ("\n")
'''



def abyss_pe_mp(pefile,mpfile):

	pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")
	mp_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_MP_Reads" + "/")
	
	with open(pefile) as fh:
		pe_sample_name_list = fh.read().splitlines()
		pe_samples=[x for x in pe_sample_name_list]
		pe_lib_string ="lib=" + "'" + ' '.join(pe_samples) + "'" + " "
		pe_left_read_name_suffix = '_R1.fastq'
		pe_right_read_name_suffix = '_R2.fastq'
		pe_left_read_name_list = [x + pe_left_read_name_suffix for x in pe_sample_name_list]
		pe_right_read_name_list = [x + pe_right_read_name_suffix for x in pe_sample_name_list]
		pe_samples=[x for x in pe_sample_name_list]
		pe_result = [sublist for sublist in zip(pe_samples, pe_left_read_name_list, pe_right_read_name_list)]
		pe_result1 =   [' '+ x +'=' + "'" + pe_cleaned_read_folder + y + " " + pe_cleaned_read_folder + z + "'" for x, y,z in pe_result]

	with open(mpfile) as fh:
		mp_sample_name_list = fh.read().splitlines()
		mp_samples=[x for x in mp_sample_name_list]
		mp_lib_string ="mp=" + "'" + ' '.join(mp_samples) + "'" + " "
		mp_left_read_name_suffix = '_R1.fastq'
		mp_right_read_name_suffix = '_R2.fastq'
		mp_left_read_name_list = [x + mp_left_read_name_suffix for x in mp_sample_name_list]
		mp_right_read_name_list = [x + mp_right_read_name_suffix for x in mp_sample_name_list]
		mp_samples=[x for x in mp_sample_name_list]
		mp_result = [sublist for sublist in zip(mp_samples, mp_left_read_name_list, mp_right_read_name_list)]
		mp_result1 =   [' '+ x +'=' + "'"+ mp_cleaned_read_folder + y + " " + mp_cleaned_read_folder + z + "'" for x, y,z in mp_result]

		pe_result2 = ' '.join(pe_result1)
		mp_result2 = ' '.join(mp_result1)
	   
		parse_string = pe_lib_string + mp_lib_string + pe_result2 + mp_result2

		return parse_string
'''
pefile = "./pe_samples.list"
mpfile = "./mp_samples.list"

abyss_input = abyss_pe_mp(pefile,mpfile)

print ("Input PE-MP Libraries for ABYSS2 Assembler")
print ("----------------------------------------------------------------------------------------")
print (abyss_input)
print ("\n")
print ("\n")
'''

def platanus_formater(samplefile):
	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'

		left_read_name_list = [ x + left_read_name_suffix for x in sample_name_list]
		right_read_name_list =[ x + right_read_name_suffix for x in sample_name_list]

		pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")

		
		result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]

		result1 = [pe_cleaned_read_folder + x + " " + pe_cleaned_read_folder + y + " " for x, y in result]
		parse_string = ''.join(result1)
		return parse_string

		parse_string = "-f " + left_reads + " " + right_reads 
		return parse_string

'''
infile = "./pe_samples.list"

platanus_peinput = platanus_formater(infile)
print ("Input Read Format for PLATANUS PE")
print ("----------------------------------------------------------------------------------------")
print (platanus_peinput)
print ("\n")
'''


def trinity_formater(samplefile):
	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'

		left_read_name_list = [ x + left_read_name_suffix for x in sample_name_list]
		right_read_name_list =[ x + right_read_name_suffix for x in sample_name_list]

		left_reads = ','.join(left_read_name_list)
		right_reads = ','.join(right_read_name_list)

		return "--left " + left_reads + " --right " + right_reads


'''
infile = "./pe_samples.list"

trinity_peinput = trinity_formater(infile)
print ("Input Read Format for Trinity PE")
print ("----------------------------------------------------------------------------------------")
print (trinity_peinput)
print ("\n")
'''


def trinity_formater(samplefile):
	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		read_name_suffix = '.fastq'
		

		read_name_list = [ x + read_name_suffix for x in sample_name_list]
		
		reads = ','.join(read_name_list)
		

		return "--single " + reads


'''
infile = "./pe_samples.list"

trinity_input = trinity_formater(infile)
print ("Input Read Format for Trinity SE")
print ("----------------------------------------------------------------------------------------")
print (trinity_input)
'''