import re, random, string, os, subprocess
from general import randomstring, try_to_remove_path

def RNAfold_call(sequence,unique_string='',delete = False,extra_parameter_set_configuration='',noPS=True): #calls RNAfold and returns stdout
	if unique_string == '':
		unique_string = randomstring(20)
	if noPS == True:
		PS = '--noPS'
	elif noPS == False:
		PS = ''
	current_dir = os.getcwd()
	sequence = translate_to_RNA(sequence)
	p = subprocess.Popen(['printf ">'+unique_string+'\n'+sequence+'" | RNAfold -p '+PS+' --noLP -d2 --bppmThreshold=1e-100'+extra_parameter_set_configuration], shell=True, stdout=subprocess.PIPE)
	p.wait()
	if delete is True:
		try_to_remove_path(current_dir+'/'+unique_string+'_dp.ps')
	return p.stdout.read()
	
def extract_dG(fold_output,delete=False): #extracts the Gibbs free energy of folds
	MFE = re.search(r'\(\ *(\-*\ *\d{1,2}\.\d{2}).*',fold_output)
	dG=MFE.group(1)
	if delete == True:
		path_re = re.search(r'>(.+)\n',fold_output)
		path_fold_output = os.getcwd()+'/'+path_re.group(1)+'_dp.ps'
		try_to_remove_path(path_fold_output)
	return float(dG)

def extract_sec_str(fold_output):
	sec_str_re = re.search(r'([\.\(\)]+)\ ',fold_output)
	sec_str = sec_str_re.group(1)
	return sec_str
	
def extract_dG_ensemble(fold_output,delete=False): #extracts the Gibbs free energy of folds
	MFE = re.search(r'\[\ *(\-*\ *\d{1,2}\.\d{2}).*',fold_output)
	if MFE:
		dG=MFE.group(1)
	else:
		dG=0
	if delete == True:
		path_re = re.search(r'>(.+)\n',fold_output)
		path_fold_output = os.getcwd()+'/'+path_re.group(1)+'_dp.ps'
		try_to_remove_path(path_fold_output)
	return float(dG)

def translate_to_RNA(sequence):#turns DNA string into RNA (upper case)
	sequence_as_RNA=str(sequence).translate(string.maketrans("Tatgcunwsmkrybdhv", "UAUGCUNWSMKRYBDHV")) #removes all bad characters
	return sequence_as_RNA
	
def extract_probabilities_from_dot_plot(path,range_of_nucleotides,delete=False): #nucleotides defined in python style: first nucleotide = 0
	assert isinstance(range_of_nucleotides,list)
	nucleotides_list = []
	probability_list = []
	end_nucleotide = range_of_nucleotides[1]+2
	#Checking if end_nucleotide isnt higher than sequence used in RNAfold:
	f = open(path,'r')
	for line in f:
		sequence_found = re.search(r'([AGUC]+)\\',line)
		if sequence_found:
			if len(sequence_found.group(1))+1<end_nucleotide:
				end_nucleotide = len(sequence_found.group(1))+1
	for nucleotide in range(range_of_nucleotides[0]+1,end_nucleotide):
		probability = float(0)
		nucleotide = str(nucleotide)
		f = open(path,'r')
		for line in f:
			result_of_search2 = re.search(r"^"+nucleotide+"\s\d+\s(\d+\.\d+)\subox",line)
			if result_of_search2:
				probability += float(result_of_search2.group(1))**2
			result_of_search = re.search(r"^\d+\s"+nucleotide+"\s(\d+\.\d+)\subox",line)
			if result_of_search:
				probability += float(result_of_search.group(1))**2
		probability_list.append(probability)
		nucleotides_list.append(int(nucleotide))
	if delete == True:
		try_to_remove_path(path)
	return (nucleotides_list,probability_list)
	
def RNAfold_constrained_call(sequence,constraint,unique_string='',extra_parameter_set_configuration='',noPS=True):
	if unique_string=='':
		unique_string = randomstring(20)
	if noPS == True:
		PS = '--noPS'
	elif noPS == False:
		PS = ''
	sequence = translate_to_RNA(sequence)
	p = subprocess.Popen(['printf ">'+unique_string+'\n'+sequence+'\n'+constraint+'" | RNAfold -C -p '+PS+' --noLP -d2 --bppmThreshold=1e-100'+extra_parameter_set_configuration], shell=True, stdout=subprocess.PIPE)
	p.wait()
	return p.stdout.read()

def extract_probabilities_from_nucleotide_list(nucleotide_list,source_nucleotide_list,source_probability_list):
	prob_list = []
	for i in range(len(source_nucleotide_list)):
		if source_nucleotide_list[i] in nucleotide_list:
			prob_list.append(source_probability_list[i])
	return(prob_list)
	
def RNAcofold_call(seq1,seq2,unique_string='',extra_parameter_set_configuration='',noPS=True): #calls RNAfold and returns stdout
	seq1 = translate_to_RNA(seq1)
	seq2 = translate_to_RNA(seq2)
	if noPS == True:
		PS = '--noPS'
	elif noPS == False:
		PS = ''
	if unique_string == '':
		unique_string=randomstring(20)
	p = subprocess.Popen(['printf ">'+unique_string+'\n'+seq1+'&'+seq2+'" | RNAcofold -p '+PS+' --noLP -d2 --bppmThreshold=1e-100'+extra_parameter_set_configuration], shell=True, stdout=subprocess.PIPE)
	p.wait()
	path_fold = os.getcwd()
	return p.stdout.read()

def RNAcofold_call_constrained(seq1,seq2,constr1,constr2,unique_string='',extra_parameter_set_configuration='',noPS=True): #calls RNAcofold and returns stdout
	seq1 = translate_to_RNA(seq1)
	seq2 = translate_to_RNA(seq2)
	if noPS == True:
		PS = '--noPS'
	elif noPS == False:
		PS = ''
	if unique_string == '':
		unique_string=randomstring(20)
	p = subprocess.Popen(['printf ">'+unique_string+'\n'+seq1+'&'+seq2+'\n'+constr1+'&'+constr2+'" | RNAcofold -C -p '+PS+' --noLP -d2 --bppmThreshold=1e-100'+extra_parameter_set_configuration], shell=True, stdout=subprocess.PIPE)
	p.wait()
	return p.stdout.read()
	 
def extract_average_probability_from_nucleotide_list(nucleotide_list,source_nucleotide_list,source_probability_list):
	prob_list = []
	for i in range(len(source_nucleotide_list)):
		if source_nucleotide_list[i] in nucleotide_list:
			prob_list.append(source_probability_list[i])
	return(sum(prob_list)/len(prob_list))
	
def find_number_of_AUG(input_sequence):
	input_sequence = translate_to_RNA(input_sequence)
	return(len([m.start() for m in re.finditer('AUG',input_sequence)]))

def extract_all_nucleotide_binding_probabilities(path_file_partition_function,RNAcofold=False,delete=False):
	path_re = re.search(r'>(.+)\n',path_file_partition_function)
	if path_re:
		path_file_partition_function = os.getcwd()+'/'+path_re.group(1)+'_dp.ps'
	if RNAcofold is True:
		pf_output_fh = open(path_file_partition_function,'r')
		full_output = pf_output_fh.read()
		cutoff_point = int(re.search(r'/cutpoint\s(\d+)\sdef',full_output).group(1))
		sequence_re = re.search(r'\/sequence { \(\\\n([AGUC]+)\\',full_output)
		sequence_full = sequence_re.group(1)
		seq1 = sequence_full[:cutoff_point-1]
		nr_seq1 = [int(x+1) for x in range(len(seq1))]
		seq2 = sequence_full[cutoff_point-1:]
		nr_seq2 = [int(x+cutoff_point) for x in range(len(seq2))]
		p_inter_seq1 = {i:float(0) for i in nr_seq1}
		p_intra_seq1 = {i:float(0) for i in nr_seq1}
		p_inter_seq2 = {i:float(0) for i in nr_seq2}
		p_intra_seq2 = {i:float(0) for i in nr_seq2}
		#[float(0)] * l_seq1
		lines_with_probabilities = re.finditer(r'(?=\n(\d+)\s(\d+)\s(\d\.\d+)\subox)',full_output)
		for match in lines_with_probabilities:
			nr1 = int(match.group(1))
			nr2 = int(match.group(2))
			prob = float(match.group(3))**2
		#	print str(nr1)+'  --  '+str(nr2)
			if nr1 in nr_seq1:
				if nr2 in nr_seq1:
			#		print 'intra'
			#		print '1 -> 1 : '+str(prob)
			#		print str(nr1)+' -> '+str(nr2)+' : '+str(prob)	
					p_intra_seq1[nr1] += prob
					p_intra_seq1[nr2] += prob
				elif nr2 in nr_seq2:
			#		print 'inter'
			#		print '1 -> 2 : '+str(prob)
			#		print str(nr1)+' -> '+str(nr2)+' : '+str(prob)
					p_inter_seq1[nr1] += prob
					p_inter_seq2[nr2] += prob
			elif nr1 in nr_seq2:
				if nr2 in nr_seq1:
			#		print 'inter'
			#		print '2 -> 1 : '+str(prob)
			#		print str(nr1)+' -> '+str(nr2)+' : '+str(prob)
					p_inter_seq2[nr1] += prob
					p_inter_seq1[nr2] += prob
				elif nr2 in nr_seq2:
			#		print 'intra'
			#		print '2 -> 2 : '+str(prob)
			#		print str(nr1)+' -> '+str(nr2)+' : '+str(prob)
					p_intra_seq2[nr1] += prob
					p_intra_seq2[nr2] += prob
	#	p_intra_seq1_list =  [p_intra_seq1[i] for i in nr_seq1]
	#	p_intra_seq2_list =  [p_intra_seq2[i] for i in nr_seq2]
	#	p_inter_seq1_list =  [p_inter_seq1[i] for i in nr_seq1]
	#	p_inter_seq2_list =  [p_inter_seq2[i] for i in nr_seq2]
		return p_intra_seq1, p_intra_seq2, p_inter_seq1, p_inter_seq2
	elif RNAcofold is False:
		pf_output_fh = open(path_file_partition_function,'r')
		full_output = pf_output_fh.read()
		sequence_re = re.search(r'\/sequence { \(\\\n([AGUC]+)\\',full_output)
		sequence_full = sequence_re.group(1)
		nr_sequence = [int(x+1) for x in range(len(sequence_full))]
		p_intra_seq = {i:float(0) for i in nr_sequence}
		lines_with_probabilities = re.finditer(r'(?=\n(\d+)\s(\d+)\s(\d\.\d+)\subox)',full_output)
		for match in lines_with_probabilities:
			nr1 = int(match.group(1))
			nr2 = int(match.group(2))
			prob = float(match.group(3))**2
			p_intra_seq[nr1] += prob
			p_intra_seq[nr2] += prob
			
		if delete is True:
			try_to_remove_path(path_file_partition_function)
		#	print str(nr1)+'  --  '+str(nr2)+'   \t   :'+str(prob)
		return p_intra_seq

def randomrna(length): #generates random strings of RNA with length given as argument
	bases = ['A','U','G','C']
	return ''.join(random.choice(bases) for i in range(length))
			
def RNAsubopt_call(sequence,n,extra_parameter_set_configuration=''): #returns subopt call
	unique_string = randomstring(20)
	sequence = translate_to_RNA(sequence)
	current_dir = os.getcwd()
	p = subprocess.Popen(['printf ">'+unique_string+'\n'+sequence+'" | RNAsubopt -p '+str(n)+extra_parameter_set_configuration], shell=True, stdout=subprocess.PIPE)
	p.wait()
#	if delete == True:
#		os.remove(current_dir+'/'+unique_string+'_dp.ps')
	return p.stdout.read()


def RNAup_call(seq1,seq2,extra_parameter_set_configuration=''): #returns RNAup call
	seq1 = translate_to_RNA(seq1)
	seq2 = translate_to_RNA(seq2)
	p = subprocess.Popen(['echo "'+seq1+'&'+seq2+'" | RNAup --no_output_file  -b'+extra_parameter_set_configuration], shell=True, stdout=subprocess.PIPE)
	p.wait()
	return p.stdout.read()
		
def extract_int_tot_seed(RNAup_fold): #extracts interaction and total term of RNAup
	m = re.search(r'\s\(([-]*\d+\.\d+)\s=\s([-]*\d+\.\d+)',RNAup_fold)
	return(float(m.group(1)),float(m.group(2)))
