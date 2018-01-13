import re, string, subprocess, random, os, itertools, os.path

def RNAfold_call(sequence,path='',unique_string=''): #calls RNAfold and returns stdout
	if unique_string == '':
		unique_string = randomstring(25)
	sequence = translate_to_RNA(sequence)
	if len(path)>0:
		p = subprocess.Popen(['cd '+path+' && printf ">'+unique_string+'\n'+sequence+'" | RNAfold -p --noPS --noLP -d2 --bppmThreshold=1e-100'], shell=True, stdout=subprocess.PIPE)
	else:
		p = subprocess.Popen(['printf ">'+unique_string+'\n'+sequence+'" | RNAfold -p --noPS --noLP -d2 --bppmThreshold=1e-100'], shell=True, stdout=subprocess.PIPE)
	p.wait()
	return p.stdout.read()
	
def RNAcofold_call(seq1,seq2,path='',unique_string='unique_string_RNAcofold'): #calls RNAcofold and returns stdout
	seq1 = translate_to_RNA(seq1)
	seq2 = translate_to_RNA(seq2)
	if len(path)>0:
		p = subprocess.Popen(['cd '+path+' && printf ">'+unique_string+'\n'+seq1+'&'+seq2+'" | RNAcofold -p --noPS --noLP -d2 --bppmThreshold=1e-100'], shell=True, stdout=subprocess.PIPE)
	else:
		p = subprocess.Popen(['printf ">'+unique_string+'\n'+seq1+'&'+seq2+'" | RNAcofold -p --noPS --noLP -d2 --bppmThreshold=1e-100'], shell=True, stdout=subprocess.PIPE)
	p.wait()
	return p.stdout.read()
	
	
def randomstring(length): # requires	import random
	return ''.join(random.choice(string.digits+string.lowercase) for i in range(length))
	
def translate_to_RNA(sequence):#turns DNA string into RNA (upper case)
	sequence_as_RNA=sequence.translate(string.maketrans("Tatgcun", "UAUGCUN")) #removes all bad characters
	return sequence_as_RNA

	
def extract_dG_ensemble(fold_output,delete=False): #extracts the Gibbs free energy of folds
	MFE = re.search(r'\[\ *(\-*\ *\d{1,2}\.\d{2}).*',fold_output)
	dG=MFE.group(1)
	if delete == True:
		path_re = re.search(r'>(.+)\n',fold_output)
		path_fold_output = os.getcwd()+'/'+path_re.group(1)+'_dp.ps'
		try_to_remove_path(path_fold_output)
	return float(dG)
	
def try_to_remove_path(path):
	try:
		os.remove(path)
	except OSError:
		pass
		
def analyze_UTR(utr,cds):
	utr=translate_to_RNA(utr)
	cds=translate_to_RNA(cds)	
	cds=cds[:50]
	dG_EFE = extract_dG_ensemble(RNAfold_call(utr+cds),delete=True)
	purineAG_in_min3 = 1 if utr[-3:-2] in ['A','G'] else -1
	U_in_min3 = 1 if utr[-3:-2] == 'U' else -1
	A_in_min1 = 1 if utr[-1:] == 'A' else -1
	AA_in_min32 = 1 if utr[-3:-1] == 'AA' else -1
	CG_in_min32 = 1 if utr[-3:-1] == 'CG' else -1
	AC_in_min21 = 1 if utr[-2:] == 'AC' else -1
	oof_uAUG = count_oof_uAUG(utr)
	GACA_kmer = pattern_finder('GACA',utr)
	GG_kmer = pattern_finder('GG',utr)
	CACC_kmer = pattern_finder('CACC',utr)
	CA_in_min76 = 1 if utr[-7:-5] == 'CA' else -1
	CC_in_min76 = 1 if utr[-7:-5] == 'CC' else -1
	return (dG_EFE,purineAG_in_min3,U_in_min3,A_in_min1,AA_in_min32,CG_in_min32,AC_in_min21,oof_uAUG,GACA_kmer,GG_kmer,CACC_kmer,CA_in_min76,CC_in_min76)
		
def pattern_finder(pattern,utr):
	matched = -1
	match_pattern = '.*'+pattern+'.*'
	if re.match(match_pattern,utr):
		matched = 1
	return(matched)

def count_oof_uAUG(utr,start_codon='AUG'):
	counter_oof_uAUG = 0
	interesting_positions = [i for i in [i+1 for i in range(len(utr))] if i % 3 != 0 and i > 2]
	for position in interesting_positions:
		if utr[-position:-(position-3)] == start_codon:
			counter_oof_uAUG+=1
	return(counter_oof_uAUG)
	
def write_tuple_to_file(path,data_to_write):
	fh = open(path,'a')
	fh.write(",".join(str(item) for item in data_to_write)+'\n')
		
identification = 'id'
unique_string = str(identification)+'_'+randomstring(10)
pre_random = 'aaaacaa'

CDS = 'ATGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAATGGTCACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTGATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTAGGTTATGGTTTGCAATGTTTTGCTAGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGAACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTTAATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAACTATAACTCTCACAATGTTTACATCACTGCTGACAAACAAAAGAATGGTATCAAAGCTAACTTCAAAATTAGACACAACATTGAAGATGGTGGTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCAGTCTTGTTACCAGACAACCATTACTTATCCTATCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGACCACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATTGATGAATTGTACAAATAA'
CDS_cutoff = CDS[:50]
# Enter path where you want to save the 'output_analysis' file
path_output = '/output_analysis.csv'
try_to_remove_path(path_output)
first_line = 1
# Enter path where you saved the 'Dvir_UTRs' file
fh = open('/Dvir_UTRs.csv', 'r')

counter = 0
for line in fh:
	counter +=1
	split_line = line.split(',')
	utr_random_part = split_line[0]
	protein_abundance = split_line[1].rstrip()
	utr = pre_random + utr_random_part
	utr = translate_to_RNA(utr)
	analyzed_UTR_tuple = analyze_UTR(utr,CDS_cutoff)
	if first_line == 1:
		write_tuple_to_file(path_output,('UTR_variant','protein_abundance','full_utr','dG_EFE','purineAG_in_min3','U_in_min3','A_in_min1','AA_in_min32','CG_in_min32','AC_in_min21','oof_uAUG','GACA_kmer','GG_kmer','CACC_kmer','CA_in_min76','CC_in_min76'))
		first_line = 0
	write_tuple_to_file(path_output,(utr_random_part,protein_abundance,utr,)+analyzed_UTR_tuple)
	print counter
