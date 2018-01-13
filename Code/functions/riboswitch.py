import subprocess, re, random, string, os, math
from RNA import extract_sec_str, RNAfold_call, translate_to_RNA, extract_dG, extract_dG_ensemble, extract_average_probability_from_nucleotide_list, extract_probabilities_from_dot_plot, RNAfold_constrained_call, extract_probabilities_from_nucleotide_list, RNAcofold_call, extract_all_nucleotide_binding_probabilities, RNAcofold_call_constrained
from general import randomstring, try_to_remove_path, from_features_to_response_of_PLS
from antisense import calculate_paired_termini

def mutation(candidate,list_degenerate_positions,dict_optional_nucleotides,list_available_nucleotides,number_of_replacement_list_biased):
	#replacement	
	n_inserted = len([x for x in dict_optional_nucleotides.values() if x == 1])
	n_place_for_inserts = len([x for x in dict_optional_nucleotides.values() if x == 0])
	list_degenerate_positions_in_candidate = [x for x in list_degenerate_positions if candidate[x] != '']
	n_replacement = len(list_degenerate_positions_in_candidate)
	p_insertion = min(float(n_place_for_inserts)/sum([float(n_inserted),float(n_replacement),float(n_place_for_inserts)]),0.1)
	p_deletion = min(float(n_inserted)/sum([float(n_inserted),float(n_replacement),float(n_place_for_inserts)]),0.1)
	possible_mutations = ['replacement']*int(round((float(1)-p_insertion-p_deletion)*100))+['deletion']*int(round(p_deletion*100))+['insertion']*int(round(p_insertion*100)) #replace
	new_candidate = candidate[:]
	type_of_mutation = random.choice(possible_mutations)
	if type_of_mutation == 'replacement':
		length_inserts = len(list_degenerate_positions_in_candidate)
	#	print number_of_replacement_list_biased
		indexes_to_replace = random.sample(list_degenerate_positions_in_candidate,min(random.choice(number_of_replacement_list_biased),n_replacement))
		for index_to_replace in indexes_to_replace:
			new_candidate[index_to_replace] = random.choice([x for x in list_available_nucleotides[index_to_replace] if x not in candidate[index_to_replace]])
	elif type_of_mutation == 'deletion':
		random_nucleotide_pos = random.choice([x for x in dict_optional_nucleotides.keys() if dict_optional_nucleotides[x] == 1])
	#	print 'picked from'
	#	print [x for x in dict_optional_nucleotides.keys() if dict_optional_nucleotides[x] == 1]
		new_candidate[random_nucleotide_pos] = ''
	#	print type_of_mutation
	elif type_of_mutation == 'insertion':
		random_nucleotide_pos = random.choice([x for x in dict_optional_nucleotides.keys() if dict_optional_nucleotides[x] == 0])
		new_candidate[random_nucleotide_pos] = random.choice(list_available_nucleotides[random_nucleotide_pos])
	#	print type_of_mutation
	for pos_nucleotide in dict_optional_nucleotides.keys():
		if new_candidate[pos_nucleotide] == '':
			dict_optional_nucleotides[pos_nucleotide] = 0
		else:
			dict_optional_nucleotides[pos_nucleotide] = 1
	return new_candidate,dict_optional_nucleotides
		
def analyze_input(input_utr):
	start_nucleotide = 0
	degenerate_dict = {"N":["A","C","G","U"],"V":["A","C","G"],"H":["A","C","U"],"D":["A","G","U"],"B":["C","G","U"],"Y":["C","U"],"R":["A","G"],"K":["G","U"],"M":["A","C"],"S":["C","G"],"W":["A","U"],"A":"A","C":"C","G":"G","U":"U"}
	list_available_nucleotides = []
	list_if_nucleotides_are_mandatory = []
	found_randomized_parts = re.finditer(r'([\(NWSMKRYBDHV\)]+)',input_utr)
	if found_randomized_parts:
		for match in found_randomized_parts:
			constrained_part = input_utr[start_nucleotide:match.start(0)]
			for i in range(len(constrained_part)):
				list_available_nucleotides.append(constrained_part[i])
				list_if_nucleotides_are_mandatory.append(1)
			start_nucleotide = match.end(0)
			degenerate_part = match.group(0)
			find_facultative_nucleotides = re.search(r'([NWSMKRYBDHV]*)\(([NWSMKRYBDHV]+)\)([NWSMKRYBDHV]*)',degenerate_part)
			if find_facultative_nucleotides:
				for i in range(1,4):
					facultative_nucleotides=find_facultative_nucleotides.group(i)
					for j in range(len(facultative_nucleotides)):
						list_available_nucleotides.append(degenerate_dict[facultative_nucleotides[j]])
						if i == 2: 
							list_if_nucleotides_are_mandatory.append(0)
						else:
							list_if_nucleotides_are_mandatory.append(1)
			else:
				for i in range(len(degenerate_part)):
					list_available_nucleotides.append(degenerate_dict[degenerate_part[i]])
					list_if_nucleotides_are_mandatory.append(1)
	for i in range(len(input_utr[start_nucleotide:])):
			list_available_nucleotides.append(input_utr[start_nucleotide:][i])
			list_if_nucleotides_are_mandatory.append(1)
	return(list_available_nucleotides,list_if_nucleotides_are_mandatory)

def generate_initial_candidate(list_available_nucleotides,list_if_nucleotides_are_mandatory):
	list_current_candidate = []
	for i in range(len(list_if_nucleotides_are_mandatory)):
		if list_if_nucleotides_are_mandatory[i] == 1:
			list_current_candidate.append(random.choice(list_available_nucleotides[i]))
		elif list_if_nucleotides_are_mandatory[i] == 0:
			list_current_candidate.append('')
	return list_current_candidate	

def evaluate_trans_riboswitch(candidate,utr,dG_utr,cds,terminator,constrained_formation,energy_difference,switching_nucleotides,print_folds=False):
	cds = cds[:50]
	sequence = ''.join(candidate)
	full_utr = utr+cds
	utr_constraint = '.'*len(full_utr)
	RNAcofold = RNAcofold_call(sequence+terminator,full_utr)
	RNAcofold_constrained = RNAcofold_call_constrained(sequence+terminator,full_utr,constrained_formation,utr_constraint)
	if print_folds==True:
		print RNAcofold
		print RNAcofold_constrained
	p_intra_seq1, p_intra_seq2, p_inter_seq1, p_inter_seq2 = extract_all_nucleotide_binding_probabilities(RNAcofold,RNAcofold = True)
	p_intra_seq1_constr, p_intra_seq2_constr, p_inter_seq1_constr, p_inter_seq2_constr = extract_all_nucleotide_binding_probabilities(RNAcofold_constrained,RNAcofold = True)
	minimum_utr = min(p_inter_seq2.keys())
	converted_switching_nucleotides = [(x+minimum_utr-1) for x in switching_nucleotides]
	average_nucleotide_switching_OFF = float(sum([p_inter_seq2[x] for x in converted_switching_nucleotides]))/float(len(converted_switching_nucleotides))
	average_nucleotide_switching_ON = float(sum([p_inter_seq2_constr[x] for x in converted_switching_nucleotides]))/float(len(converted_switching_nucleotides))
	average_nucleotide_UTR_intra_ON = float(sum([p_intra_seq2_constr[x] for x in converted_switching_nucleotides]))/float(len(converted_switching_nucleotides))
	maximum_utr = max(p_inter_seq2.keys())
	p_cds_mfe = [p_inter_seq2[x] for x in range(maximum_utr-50+1,maximum_utr+1)]
	p_cds_cnstr = [p_inter_seq2_constr[x] for x in range(maximum_utr-50+1,maximum_utr+1)]
	average_nucleotide_CDS_coverage_OFF = float(sum(p_cds_mfe)/len(p_cds_mfe))
	average_nucleotide_CDS_coverage_ON = float(sum(p_cds_cnstr)/len(p_cds_cnstr))
	score_CDS_pairing = average_nucleotide_CDS_coverage_OFF+average_nucleotide_CDS_coverage_ON
	score_pairing = average_nucleotide_switching_OFF-average_nucleotide_switching_ON
	dG_AB_ON_ensemble = extract_dG_ensemble(RNAcofold_constrained,delete=True)
	dG_AB_OFF_ensemble = extract_dG_ensemble(RNAcofold,delete=True)
	
	dG_ON_ensemble = extract_dG_ensemble(RNAfold_constrained_call(sequence+terminator,constrained_formation),delete=True)
	dG_OFF_ensemble = extract_dG_ensemble(RNAfold_call(sequence+terminator),delete=True)
	score_energy_diff_aptamer_formation = dG_ON_ensemble-dG_OFF_ensemble
	ddG_ON = dG_AB_ON_ensemble-dG_ON_ensemble-dG_utr #should be maximized
	ddG_OFF = dG_AB_OFF_ensemble-dG_OFF_ensemble-dG_utr
	score_energy_formation = abs(ddG_OFF-ddG_ON-energy_difference)
	paired_termini = calculate_paired_termini(sequence,utr)/(len(sequence)/2)
	return score_energy_formation,score_pairing,average_nucleotide_switching_OFF,average_nucleotide_switching_ON,average_nucleotide_UTR_intra_ON,dG_OFF_ensemble,dG_ON_ensemble,paired_termini,score_CDS_pairing,score_energy_diff_aptamer_formation

def determine_switching_nucleotides(initial_constraint_utr):
	switching_nucleotides = []
	for nucleotide in range(1,len(initial_constraint_utr)+1):#R style vector
		if initial_constraint_utr[nucleotide-1] == '|':
			switching_nucleotides.append(nucleotide)
	return switching_nucleotides
	
def generate_constrained_form(candidate_list,aptamer_or_expression_part_list,aptamer_fold,terminator): #generates constrained string for RNA(co)fold input based on list with candidate (candidate_list)
	constrained_formation = ''
	aptamer_fold_constrained = aptamer_fold.translate(string.maketrans('.','x'))
	aptamer_constrained_list = ['aptamer']*len(aptamer_fold_constrained)
	for i in range(len(candidate_list)):
		if aptamer_or_expression_part_list[i] == 'expression' and candidate_list[i] in ['A','G','U','C']:
			constrained_formation += '.'
		elif set(aptamer_or_expression_part_list[i:i+len(aptamer_fold_constrained)]) == set(aptamer_constrained_list):
			constrained_formation += aptamer_fold_constrained
	constrained_formation += '.'*len(terminator)
	return constrained_formation

def dG_spacing_Salis(RBS_part_UTR,rRNA_16S3prime = 'AUUCCUCCA'):
#	print RBS_part_UTR
#	subopt -multi -T 37 Desktop/test_subopt
	path = '/home/gpeters/RBS_analysis_Salis.subopt'
	fh = open(path+'.in','w')
	fh.write("2\n"+translate_to_RNA(RBS_part_UTR)+"\n"+translate_to_RNA(rRNA_16S3prime)+"\n1 2\n2")
	fh.flush()
	p = subprocess.Popen(['subopt -multi -T 37 -material rna1999 '+path], shell=True, stdout=subprocess.PIPE)
	p.wait()
	fh_output = open(path+'.subopt')
	full_output = fh_output.read()
	split_output = full_output.split('% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n\n% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %')
	length_rRNA = len(rRNA_16S3prime)
	length_utr = len(RBS_part_UTR)
	length_re = length_utr+length_rRNA+1
	result_list = []
	dG_spacing_list = []
	dG_16SrRNA_list = []
	for i in range(len(split_output)):
	#	print split_output[i]
	#	print str(length_utr+length_rRNA)
		dG_search = re.search(r''+str(length_utr+length_rRNA)+'\n(\-*\d+\.\d+)\n',split_output[i])
		dG_16SrRNA = float(dG_search.group(1))
		binding_search = re.finditer(r'(\d+)\t(\d+)\n',split_output[i])
		n1 = 0
		for match in binding_search:
			bound1 = int(match.group(1))
			bound2 = int(match.group(2))
			if bound2 > length_utr:
				n1_candidate = bound2-length_utr
			#	print n1_candidate
				if n1_candidate>n1:
					n1 = n1_candidate
					n2 = bound1
#		print 'n1 = '+str(n1)+' n2 = '+str(n2)
		s = length_utr+1-n1-n2
		c1 = float(0.048)
		c2 = float(0.24)
		if s>5:
			dG_spacing = float(c1*float(s-5))**2+float(c2*(s-5))
		else:
			dG_spacing = c1/float(1+math.exp(c2*(s-5+2)))**3
		result_list.append(dG_16SrRNA+dG_spacing)
		dG_spacing_list.append(dG_spacing)
		dG_16SrRNA_list.append(dG_16SrRNA)
	index_min = result_list.index(min(result_list))
#	print 'index = '+str(index_min)
	dG_spacing = dG_spacing_list[index_min]
	return dG_spacing

def generate_dict_optional_nucleotides(list_if_nucleotides_are_mandatory):
	dict_optional_nucleotides = dict()
	for i in range(len(list_if_nucleotides_are_mandatory)):
		if list_if_nucleotides_are_mandatory[i] == 0:
			dict_optional_nucleotides[i] = 0
	return dict_optional_nucleotides

def list_degenerate_positions(list_available_nucleotides):
	list_degenerate_positions = []
	for i in range(len(list_available_nucleotides)):
		if len(list_available_nucleotides[i]) > 1:
			list_degenerate_positions.append(i)
	return list_degenerate_positions

def generate_aptamer_sequence__constrained_list_and_constrained_sec_str_aptamer(rs_constraint_original,aptamer):
	rs_constraint = rs_constraint_original.translate(string.maketrans("Tatgcunx", "UAUGCUNX"))
	aptamer_or_expression_part_list = [] #fill with expression for degenerated part (parts other than aptamer part), fixed for fixed part (aptamer for instance)
	for i in range(len(rs_constraint)):
		if rs_constraint[i] in ['N','W','S','M','K','R','Y,','B','D','H','V','A','G','U','C']:
			aptamer_or_expression_part_list.append('expression')
		elif rs_constraint[i] == 'X':
			aptamer_or_expression_part_list.append('aptamer')
	rs_constraint_re = re.search(r'([\(NWSMKRYBDHVAGCU\)]*)(X+)([\(NWSMKRYBDHVAGUC\)]*)',rs_constraint)
	if rs_constraint_re:
		preseq = rs_constraint_re.group(1)
		postseq = rs_constraint_re.group(3)
		aptamer_place = rs_constraint_re.group(2)
		rs_with_aptamer_sequence_constraint = preseq+aptamer+postseq
		if len(aptamer_place) != len(aptamer):
			raise ValueError('Length aptamer sequence and length riboswitch constraint not the same.')
	aptamer_fold = extract_sec_str(RNAfold_call(aptamer,delete=True))
	return(rs_with_aptamer_sequence_constraint,aptamer_fold,aptamer_or_expression_part_list)

def make_dict_replacement_options(list_degenerate_positions,list_available_nucleotides):
	dict_replacement_options = dict() #dict with replacements for each possition (contains another dict with possible replacements with current nucleotide as key)
	for i in list_degenerate_positions:
		dict_possible_replacement = dict()
		for possible_nucleotide in list_available_nucleotides[i]:
			list_possible_replacement = [x for x in list_available_nucleotides[i] if x != possible_nucleotide]
			dict_possible_replacement[possible_nucleotide] = list_possible_replacement
		dict_replacement_options[i] = dict_possible_replacement
	return dict_replacement_options


def create_list_with_number_of_replacement(list_degenerate_positions,list_if_nucleotides_are_mandatory):
	number_of_replacement_list_bias_to1 = []
	list_temp = []
	for x in range(len(list_degenerate_positions)):
		if list_if_nucleotides_are_mandatory[x] == 1:
			list_temp.append(1)
	for i in range(0,int(len(list_temp)/1.5)):
		i = i+1
		number_of_replacement_list_bias_to1 = number_of_replacement_list_bias_to1 + range(1,i)
	return number_of_replacement_list_bias_to1

def predict_AR_through_PLS(rs_seq,names_features,theo_apt,apt_constrained,cds,coefficient_list,scaling_list,intercept,verbose=False):
	features_rs_dict = analyze_features_RS_PLS_theo_mkate2(rs_seq,names_features,theo_apt,apt_constrained,cds,verbose=verbose)
	features_rs_list = [features_rs_dict[key] for key in names_features]
	if verbose == True:
		print features_rs_list
		print features_rs_dict
	predicted_AR_log = from_features_to_response_of_PLS(features_rs_list,names_features,coefficient_list,scaling_list,intercept)
	return predicted_AR_log


def analyze_features_RS_PLS_theo_mkate2(RS_seq,list_of_features_to_analyze,theo_apt,apt_constrained,CDS = 'AUGGUUAGCGAG',unique_string='',verbose=False):
	features_output = dict()
	if unique_string == '':
		unique_string = randomstring(20)
	CDS = CDS[:12]
	RS_seq = translate_to_RNA(RS_seq)+CDS
	begin_apt = RS_seq.find(theo_apt)
	end_apt = begin_apt+len(theo_apt)
	sec_str_apt = '.'*begin_apt+apt_constrained+'.'*(len(RS_seq)-end_apt)
	fold_mfe = RNAfold_call(RS_seq)#calls RNAfold to get MFE
	sec_str_mfe = extract_sec_str(fold_mfe)
	if verbose == True:
		print 'MFE str. = '+sec_str_mfe
	if 'dG_mfe' in list_of_features_to_analyze:
		dG_mfe = extract_dG(fold_mfe)
		if verbose == True:
			print dG_mfe
		features_output['dG_mfe'] = dG_mfe
	path_re = re.search(r'>(.+)\n',fold_mfe)
	if 'dG_mfe_ensemble' in list_of_features_to_analyze:
		dG_mfe_ensemble = extract_dG_ensemble(fold_mfe)
		if verbose == True:
			print 'dG MFE = '+str(dG_mfe_ensemble)
		features_output['dG_mfe_ensemble'] = dG_mfe_ensemble
	path_fold_output = os.getcwd()+'/'+path_re.group(1)+'_dp.ps'
	nucleotides_list_mfe,probability_list_mfe = extract_probabilities_from_dot_plot(path_fold_output,[0,len(sec_str_mfe)-1],delete=True)	
	fold_constrained=RNAfold_constrained_call(RS_seq,sec_str_apt)
	path_re = re.search(r'>(.+)\n',fold_constrained)
	path_fold_output_constrained = os.getcwd()+'/'+path_re.group(1)+'_dp.ps'
	sec_str_constrained = extract_sec_str(fold_constrained)
	if verbose == True:
		print 'CST str. = '+sec_str_constrained
	if 'dG_constrained' in list_of_features_to_analyze:
		dG_constrained = extract_dG(fold_constrained)
		if verbose == True:
			print 'dG_CST = '+str(dG_constrained)
		features_output['dG_constrained'] = dG_constrained
	nucleotides_list_constrained,probability_list_constrained = extract_probabilities_from_dot_plot(path_fold_output_constrained,[0,len(sec_str_constrained)-1])
	start_CDS = RS_seq.find(translate_to_RNA(CDS))
	start_codon = [start_CDS+1,start_CDS+2,start_CDS+3]
	block_1 = [m-3 for m in start_codon]
	block_2 = [m-6 for m in start_codon]
	block_3 = [m-9 for m in start_codon]
	block_4 = [m-12 for m in start_codon]
	block_5 = [m-15 for m in start_codon]
	all_nucleotides = block_5+block_4+block_3+block_2+block_1+start_codon
	probability_list_diff = [j-i for i,j in zip(probability_list_constrained,probability_list_mfe)]
	if 'p_AUG' in list_of_features_to_analyze:
		p_AUG = extract_average_probability_from_nucleotide_list(start_codon,nucleotides_list_constrained,probability_list_diff)
		features_output['p_AUG'] = p_AUG
	if 'p_block1' in list_of_features_to_analyze:
		p_block1 = extract_average_probability_from_nucleotide_list(block_1,nucleotides_list_constrained,probability_list_diff)
		features_output['p_block1'] = p_block1
	if 'p_block2' in list_of_features_to_analyze:
		p_block2 = extract_average_probability_from_nucleotide_list(block_2,nucleotides_list_constrained,probability_list_diff)
		features_output['p_block2'] = p_block2
	if 'p_block3' in list_of_features_to_analyze:
		p_block3 = extract_average_probability_from_nucleotide_list(block_3,nucleotides_list_constrained,probability_list_diff)
		features_output['p_block3'] = p_block3
	if 'p_block4' in list_of_features_to_analyze:
		p_block4 = extract_average_probability_from_nucleotide_list(block_4,nucleotides_list_constrained,probability_list_diff)
		features_output['p_block4'] = p_block4
	if 'p_block5' in list_of_features_to_analyze:
		p_block5 = extract_average_probability_from_nucleotide_list(block_5,nucleotides_list_constrained,probability_list_diff)
		features_output['p_block5'] = p_block5
	#EXTRA ENSEMBLE ENERGIES
	all_nucleotides_CDS = range(max(all_nucleotides)+1,max(all_nucleotides)+10)
	all_probs_CDS = extract_probabilities_from_nucleotide_list(all_nucleotides_CDS,nucleotides_list_constrained,probability_list_diff)
	all_probs = extract_probabilities_from_nucleotide_list(all_nucleotides,nucleotides_list_constrained,probability_list_diff)
	all_probs = list(reversed(all_probs)) #REVERSE TO make P0 first value

	probability_names = ['p'+str(i) for i in range(0,18,1)]
	probability_names_CDS = ['p_CDS'+str(i) for i in range(1,10)]
	for i in range(len(probability_names)):
		features_output[probability_names[i]] = all_probs[i]
	for i in range(len(probability_names_CDS)):
		features_output[probability_names_CDS[i]] = all_probs_CDS[i]
	if 'dG_constrained_ensemble' in list_of_features_to_analyze:
		dG_constrained_ensemble = extract_dG_ensemble(fold_constrained,delete=True)
		if verbose == True:
			print 'dG_CST = '+str(dG_constrained_ensemble)
		features_output['dG_constrained_ensemble'] = dG_constrained_ensemble
	utr_rs_only = RS_seq[:start_CDS]
	utr_rs_only_short = RS_seq[end_apt:start_CDS]
	rRNA_16S = 'AUUCCUCCA'
	if 'dG_16S_rRNA' or 'dG_16S_rRNA_ensemble' in list_of_features_to_analyze:
		fold_16S_rRNA = RNAcofold_call(utr_rs_only,rRNA_16S)
	if 'dG_16S_rRNA' in list_of_features_to_analyze:
		if 'dG_16S_rRNA_ensemble' in list_of_features_to_analyze:
			dG_16S_rRNA = extract_dG(fold_16S_rRNA)
		else:
			dG_16S_rRNA = extract_dG(fold_16S_rRNA,delete=True)
		features_output['dG_16S_rRNA'] = dG_16S_rRNA
	if 'dG_16S_rRNA_ensemble' in list_of_features_to_analyze:
		dG_16S_rRNA_ensemble = extract_dG_ensemble(fold_16S_rRNA,delete=True)
		features_output['dG_16S_rRNA_ensemble'] = dG_16S_rRNA_ensemble
	if 'dG_16S_rRNA_constrained' or 'dG_16S_rRNA_constrained_ensemble' in list_of_features_to_analyze:
		fold_16S_rRNA_constrained = RNAcofold_call_constrained(utr_rs_only,rRNA_16S,sec_str_apt[:start_CDS],'.'*len(rRNA_16S))
	if 'dG_16S_rRNA_constrained' in list_of_features_to_analyze:
		if 'dG_16S_rRNA_constrained_ensemble' not in list_of_features_to_analyze:
			dG_16S_rRNA_constrained = extract_dG(fold_16S_rRNA_constrained,delete=True)
		else:
			dG_16S_rRNA_constrained = extract_dG(fold_16S_rRNA_constrained,delete=False)
		features_output['dG_16S_rRNA_constrained'] = dG_16S_rRNA_constrained
	#fold_16S_rRNA_constrained_spacing = RNAcofold_call_constrained(utr_rs_only_short,rRNA_16S,sec_str_apt[end_apt:start_CDS],'.'*len(rRNA_16S))
	if 'dG_16S_rRNA_constrained_ensemble' in list_of_features_to_analyze:
		dG_16S_rRNA_constrained_ensemble = extract_dG_ensemble(fold_16S_rRNA_constrained,delete=True)
		features_output['dG_16S_rRNA_constrained_ensemble'] = dG_16S_rRNA_constrained_ensemble
	if 'dG_spacing_off' in list_of_features_to_analyze:
		dG_spacing_off =  dG_spacing_Salis(utr_rs_only)
		features_output['dG_spacing_off'] = dG_spacing_off
	if 'dG_spacing_on' in list_of_features_to_analyze:
		dG_spacing_on =  dG_spacing_Salis(utr_rs_only_short)
		features_output['dG_spacing_on'] = dG_spacing_on
	if 'dG_spacing_on' in list_of_features_to_analyze:
		dG_spacing_on =  dG_spacing_Salis(utr_rs_only_short)
		features_output['dG_spacing_on'] = dG_spacing_on
		
	if 'ddG_16SrRNA_ensemble' in list_of_features_to_analyze:
		ddG_16SrRNA_ensemble = features_output['dG_16S_rRNA_ensemble']-features_output['dG_16S_rRNA_constrained_ensemble']
		features_output['ddG_16SrRNA_ensemble'] = ddG_16SrRNA_ensemble 
	if 'ddG' in list_of_features_to_analyze:
		features_output['ddG'] = features_output['dG_mfe_ensemble']-features_output['dG_constrained_ensemble']
	
#	print dG_16S_rRNA, dG_16S_rRNA_ensemble, dG_16S_rRNA_constrained, dG_16S_rRNA_constrained_ensemble
#	print features_output
	return features_output

def check_interference_aptamer_part(full_sequence,aptamer_fold):
	nucleotides_before_aptamer,nucleotides_until_end_aptamer=get_nucleotides_to_end_and_beginning_of_aptamer(splitted_sequence)
	transcriptional_steps = get_transcriptional_steps(input_sequence,nucleotides_until_end_aptamer)
	interference_aptamer = False
	for step in transcriptional_steps:
		transcriptional_step = full_sequence[:step]
		if (extract_sec_str(RNAfold_call(transcriptional_step))[nucleotides_before_aptamer:nucleotides_until_end_aptamer] != aptamer_fold):
		 #	print extract_sec_str(RNAfold_call(transcriptional_step))
			interference_aptamer = True
	return interference_aptamer
	

