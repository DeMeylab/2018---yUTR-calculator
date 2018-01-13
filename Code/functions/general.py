import subprocess, re, random, string, os, math
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def append_to_fasta(path,sequence,name="",description="",identification=""):
	record = SeqRecord(Seq(sequence,
					   generic_rna),
					   id=str(identification), name=str(name),
					   description=str(description))
	output_handle = open(path, "a")
	SeqIO.write(record, output_handle, "fasta")
	output_handle.close()
	
def len_file(path): #calculates the number of lines in a file
	wc=subprocess.check_output(['wc', '-l', path])
	m = re.match('(\d+)\ ',wc) #matches number of lines in wc
	if m:
		n_lines = m.group(1) #extracts sequence
	return int(n_lines)
	
def try_to_remove_path(path): #remove file (path)
	try:
		os.remove(path)
	except OSError:
		pass

def from_features_to_response_of_PLS(features_list,names_features_list,coef_list,scales_list,intercept):
	for i in range(len(features_list)):
		
		features_list[i] = features_list[i]/scales_list[i]
	response = np.dot(features_list,coef_list)+intercept
	return response

def coef_scales_from_dict_to_list(coef_dict,scales_dict,names_features):
	coef_list = []
	scales_list = []
	for name in names_features:
		coef_list.append(coef_dict[name])
		scales_list.append(scales_dict[name])
	intercept = coef_dict['intercept']
	return coef_list,scales_list,intercept
		
def randomstring(length): # requires	import random
	return ''.join(random.choice(string.digits+string.lowercase) for i in range(length))

def extract_scale_and_coefficients(coef_csvfile,scales_csvfile):
	coef_fh = open(coef_csvfile,'r')
	scales_fh = open(scales_csvfile,'r')
	coef_dict = dict()
	scales_dict = dict()
	next(coef_fh)
	for coef_line in coef_fh:
		key = coef_line.rstrip().split(',')[2]
		value = coef_line.rstrip().split(',')[1]
		if key != '':
			coef_dict[key] = np.float128(value)
	next(scales_fh)
	for scales_line in scales_fh:
		key = scales_line.rstrip().split(',')[2]
		value = scales_line.rstrip().split(',')[1]
		if key != '':
			scales_dict[key] = np.float128(value)
	return coef_dict,scales_dict

def find_longest_repeat(string):
	previous = ''
	counter = 1
	counter_max = 0
	for letter in string:
		if letter == previous:
			counter+=1
		elif counter >=1:
			if counter > counter_max:
				counter_max = counter
			counter = 1
		else:
			counter = 1
		previous = letter
	return(counter_max)


def metropolis_algorithm(score,mutation_score,temperature,maximize=False):#metropolis 
	if maximize == False:
		if mutation_score<score:
			return(True)
		else:
			dE = float(mutation_score - score)
			p = min(1,float(math.exp(-dE/temperature)))
			randomfloat=random.uniform(0,1)
			if (randomfloat <= p):
				return(True)
			else:
				return(False)
	elif maximize == True:
		if mutation_score>score:
			return(True)
		else:
			dE = float(score - mutation_score)
			p = min(1,float(math.exp(-dE/temperature)))
			randomfloat=random.uniform(0,1)
			if (randomfloat <= p):
				return(True)
			else:
				return(False)

def create_path_if_needed(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
