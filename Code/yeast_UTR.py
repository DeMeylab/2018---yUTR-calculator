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
	
def RNAcofold_call(seq1,seq2,path='',unique_string='unique_string_RNAcofold'): #calls RNAfold and returns stdout
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
		
def find_other_features(utr,CDS):
	CDS_cutoff = CDS[:50]
	utr = translate_to_RNA(utr)
	len_utr = float(len(utr))
	nucleotides = ['A','G','C','U']
	other_features = []
	names_other_features = []
	#probability nucleotides
	for nucleotide in nucleotides:
		names_other_features.append('p'+nucleotide)
		other_features.append(float(utr.count(nucleotide))/len_utr)
		
	#probability R W M	
	degenerate_dict_partial = {'W':['A','U'],'R':['A','G'],'M':['A','C']}
	for key in degenerate_dict_partial.keys():
		names_other_features.append('p'+key)
		other_features.append(sum([float(utr.count(x)) for x in degenerate_dict_partial[key]])/len_utr)
	#dinucleotide frequency
	for combi in itertools.product(nucleotides,repeat=2):
		product_dinucleotides=''.join(combi)
		if product_dinucleotides != 'GG':
			names_other_features.append('p'+product_dinucleotides)
			other_features.append(float(utr.count(product_dinucleotides))/len_utr)
#trinucleotide frequency
#	for combi in itertools.product(nucleotides,repeat=3):
#		product_trinucleotides=''.join(combi)
#	print product_trinucleotides
#		names_other_features.append('p'+product_trinucleotides)
#		other_features.append(float(utr.count(product_trinucleotides))/len_utr)
	#4-mer nucleotide frequency
#	for combi in itertools.product(nucleotides,repeat=4):
#		product_4nucleotides=''.join(combi)
#	print product_trinucleotides
#		if product_4nucleotides not in ['CACC','GACA']:
#			names_other_features.append('p'+product_4nucleotides)
#			other_features.append(float(utr.count(product_4nucleotides))/len_utr)
			
	#inframe uAUG
	names_other_features.append('ifuAUG')
	other_features.append(count_if_uAUG(utr))
	#other in and outframe startcodon mutants
	for start_codon in ['UUG','CUG','GUG','AAG','ACG','AGG','AUA','AUU','AUC']:
		names_other_features.append('ifu'+start_codon)
		other_features.append(count_if_uAUG(utr,start_codon))
		names_other_features.append('oofu'+start_codon)
		other_features.append(count_oof_uAUG(utr,start_codon))
	#Panek sticky region affinity
##	sticky_regions_18S_rRNA = [translate_to_RNA(x) for x in ['AAATCAGTTATCGTTTATTTGATAGTTCCTTTA','TTTGGAAGAGATGTATTTATTAGATAAAAAATCAATGTCTTCGG','TCATAATAACTTTTCGAATCG','TTGGCCGGTCCGATTTTTTCGTGTACTGGAT','ACCAGGACTTTTACTTTGAAAAAATTAGAGTGTTCAAA','CGTATTGCTCGAATATATTAGCAT','CGAATATATTAGCATGGAATA','GGACGTTTGGTTCTATTTTGTTGGTTTCTAGGAC','GTGAAATTCTTGGATTTATTGAAGAC','GATCGGGTGGTGTTTTTTTAATGACCCACTCG']]
##	sticky_nr = 1
##	for sticky_region in sticky_regions_18S_rRNA:
##		fold_18S_rRNA = RNAcofold_call(utr+CDS_cutoff,sticky_region)
##		dG_18S_rRNA_ensemble = extract_dG_ensemble(fold_18S_rRNA)
##		
##		names_other_features.append('sticky'+str(sticky_nr))
##		other_features.append(dG_18S_rRNA_ensemble)
##		sticky_nr +=1
	#single on specific position
	for position in [x+1 for x in range(int(len_utr))]:
		for nucleotide in nucleotides:
			if position == 1 and nucleotide == 'A':
				pass
			if position == 3 and nucleotide == 'U':
				pass
			elif position == 1:
				nucleotide_this_position = 1 if utr[-1:] == nucleotide else -1
				names_other_features.append('-'+str(position)+nucleotide)
				other_features.append(nucleotide_this_position)
			else:
				nucleotide_this_position = 1 if utr[-position:-(position-1)] == nucleotide else -1
				names_other_features.append('-'+str(position)+nucleotide)
				other_features.append(nucleotide_this_position)
	#dinucleotide on specific position
	for position in [x+1 for x in range(int(len_utr)-1)]:
		for combi in itertools.product(nucleotides,repeat=2):
			product_dinucleotides=''.join(combi)
			if position == 2 and product_dinucleotides == 'AA':
				pass		
			elif position == 2 and product_dinucleotides == 'CG':
				pass
			elif position == 1 and product_dinucleotides == 'AC':
				pass
			elif position == 6 and product_dinucleotides == 'CA':
				pass
			elif position == 6 and product_dinucleotides == 'CC':
				pass
			elif position == 1:
				dinucleotide_this_position = 1 if utr[-2:] == product_dinucleotides else -1
				other_features.append(dinucleotide_this_position)
				names_other_features.append('-'+str(position)+product_dinucleotides)
			else:
				dinucleotide_this_position = 1 if utr[-(position+1):-(position-1)] == product_dinucleotides else -1
				other_features.append(dinucleotide_this_position)
				names_other_features.append('-'+str(position)+product_dinucleotides)
#degenerate on specific position 'W':['A','U'],'R':['A','G'],'M':['A','C']
	purineAG_in_min3 = 1 if utr[-3:-2] in ['A','G'] else -1
	for position in [x+1 for x in range(int(len_utr))]:
		for key in degenerate_dict_partial.keys():
			if position == 2 and key == 'R':
				pass
			elif position == 1:
				nucleotide_this_position = 1 if utr[-1:] in degenerate_dict_partial[key] else -1
				names_other_features.append('-'+str(position)+key)
				other_features.append(nucleotide_this_position)
			else:
				nucleotide_this_position = 1 if utr[-position:-(position-1)] in degenerate_dict_partial[key] else -1
				names_other_features.append('-'+str(position)+key)
				other_features.append(nucleotide_this_position)
		
	return([tuple(names_other_features),tuple(other_features)])
#UUG, CUG, GUG, AAG, ACG, AGG, AUA, AUU, and AUC
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
	
def count_if_uAUG(utr,start_codon='AUG'):
	counter_if_uAUG = 0
	interesting_positions = [i for i in [i+1 for i in range(len(utr))] if i % 3 == 0 and i > 2]
	for position in interesting_positions:
		if position == 3:
			if utr[-position:] == start_codon:
				counter_if_uAUG+=1
		else:
			if utr[-position:-(position-3)] == 'AUG':
				counter_if_uAUG+=1
	return(counter_if_uAUG)
	
def write_tuple_to_file(path,data_to_write):
	fh = open(path,'a')
	fh.write(",".join(str(item) for item in data_to_write)+'\n')
	
	
	
	
	
	
	
identification = 'id'
unique_string = str(identification)+'_'+randomstring(10)
pre_random = 'aaaacaa'

CDS = 'ATGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAATGGTCACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTGATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTAGGTTATGGTTTGCAATGTTTTGCTAGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGAACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTTAATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAACTATAACTCTCACAATGTTTACATCACTGCTGACAAACAAAAGAATGGTATCAAAGCTAACTTCAAAATTAGACACAACATTGAAGATGGTGGTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCAGTCTTGTTACCAGACAACCATTACTTATCCTATCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGACCACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATTGATGAATTGTACAAATAA'
CDS_cutoff = CDS[:50]
path_output = '/home/gpeters/Dropbox/Code/YeastUTR/output_analysis.csv'
try_to_remove_path(path_output)
#already_done = []
first_line = 1
#if os.path.isfile(path_output):
#	fh_output = open(path_output,'r')
#	first_line = 0
#	for line in fh_output:
#		already_done.append(line.split(',')[0])
fh = open('/home/gpeters/Dropbox/Code/YeastUTR/Dvir_UTRs.csv', 'r')

counter = 0
for line in fh:
	counter +=1
	split_line = line.split(',')
	utr_random_part = split_line[0]
#	if utr_random_part not in already_done:
	protein_abundance = split_line[1].rstrip()
	utr = pre_random+utr_random_part
	utr=translate_to_RNA(utr)
	analyzed_UTR_tuple = analyze_UTR(utr,CDS_cutoff)
		#print utr
	other_features = find_other_features(utr,CDS_cutoff)
	if first_line == 1:
		write_tuple_to_file(path_output,('UTR_variant','protein_abundance','full_utr','dG_EFE','purineAG_in_min3','U_in_min3','A_in_min1','AA_in_min32','CG_in_min32','AC_in_min21','oof_uAUG','GACA_kmer','GG_kmer','CACC_kmer','CA_in_min76','CC_in_min76'))#+other_features[0])
		first_line = 0
	write_tuple_to_file(path_output,(utr_random_part,protein_abundance,utr,)+analyzed_UTR_tuple)#+other_features[1])
	print counter
