from RNA import extract_all_nucleotide_binding_probabilities, RNAup_call, extract_int_tot_seed, RNAsubopt_call, translate_to_RNA, RNAcofold_call, RNAfold_call, extract_dG_ensemble
from general import randomstring, try_to_remove_path
import re, os, time, random


def mutation(candidate,list_degenerate_positions,dict_optional_nucleotides,list_available_nucleotides,number_of_replacement_list_biased):
	#replacement	
	n_inserted = len([x for x in dict_optional_nucleotides.values() if x == 1])
	n_place_for_inserts = len([x for x in dict_optional_nucleotides.values() if x == 0])
	list_degenerate_positions_in_candidate = [x for x in list_degenerate_positions if candidate[x] != '']
	n_replacement = len(list_degenerate_positions_in_candidate)
	p_insertion = min(float(n_place_for_inserts)/sum([float(n_inserted),float(n_replacement),float(n_place_for_inserts)]),0.2)#anders dan in riboswitch zoektocht via PLS
	p_deletion = min(float(n_inserted)/sum([float(n_inserted),float(n_replacement),float(n_place_for_inserts)]),0.2)#anders dan in riboswitch zoektocht via PLS
	possible_mutations = ['replacement']*int(round((float(1)-p_insertion-p_deletion)*100))+['deletion']*int(round(p_deletion*100))+['insertion']*int(round(p_insertion*100)) #replace
	new_candidate = candidate[:]
	type_of_mutation = random.choice(possible_mutations)
	if type_of_mutation == 'replacement':
		length_inserts = len(list_degenerate_positions_in_candidate)
		indexes_to_replace = random.sample(list_degenerate_positions_in_candidate,min(random.choice(number_of_replacement_list_biased),n_replacement))
		for index_to_replace in indexes_to_replace:
			new_candidate[index_to_replace] = random.choice([x for x in list_available_nucleotides[index_to_replace] if x not in candidate[index_to_replace]])
	elif type_of_mutation == 'deletion':
		random_nucleotide_pos = random.choice([x for x in dict_optional_nucleotides.keys() if dict_optional_nucleotides[x] == 1])
		new_candidate[random_nucleotide_pos] = ''
	elif type_of_mutation == 'insertion':
		random_nucleotide_pos = random.choice([x for x in dict_optional_nucleotides.keys() if dict_optional_nucleotides[x] == 0])
		new_candidate[random_nucleotide_pos] = random.choice(list_available_nucleotides[random_nucleotide_pos])
	for pos_nucleotide in dict_optional_nucleotides.keys():
		if new_candidate[pos_nucleotide] == '':
			dict_optional_nucleotides[pos_nucleotide] = 0
		else:
			dict_optional_nucleotides[pos_nucleotide] = 1
	return new_candidate,dict_optional_nucleotides
		
def AS_UTR_analysis(AS,UTR,p_UTR,CDS,real_RBS,unique_string=''):
	#print p_UTR
	CDS = translate_to_RNA(CDS)[:50]
	if unique_string == '':
		unique_string=randomstring(20)
	dG_UTR_AS = extract_dG_ensemble(RNAcofold_call(UTR+CDS,AS,unique_string))
	dp_file_path = os.getcwd()+'/'+unique_string+'_dp.ps'
	p_intra_seq1, p_intra_seq2, p_inter_seq1, p_inter_seq2 = extract_all_nucleotide_binding_probabilities(dp_file_path,RNAcofold = True)
	#print p_intra_seq1
	#print p_intra_seq2
	real_RBS5 = real_RBS[3:-3]
	p_RBS11 = []
	p_RBS5  = []
	for nucleotide in real_RBS:
		p_RBS11.append(p_inter_seq1[nucleotide])
		if nucleotide in real_RBS5:
			p_RBS5.append(p_inter_seq1[nucleotide])
	rbs_coverage_5 = sum(p_RBS5)/len(p_RBS5)
	rbs_coverage_11 = sum(p_RBS11)/len(p_RBS11)
	p_AU = []
	for nucleotide in p_UTR.keys():
	#	print str(nucleotide)+'  pUTR = '+str(p_UTR[nucleotide])+'  pAS_UTR = '+str(p_intra_seq1[nucleotide])
		p_AU.append((1-p_UTR[nucleotide])*p_inter_seq1[nucleotide])
	p_U_avail = sum(p_AU)/sum(p_UTR.values())
	try_to_remove_path(dp_file_path)
#	print 'einde'
	return(dG_UTR_AS,p_U_avail,rbs_coverage_5,rbs_coverage_11)

def AS_analysis(features_list,AS,terminator,UTR,dG_UTR,p_UTR,real_RBS,seeds_collection,dict_utr_probs,CDS='',extra_parameter_set_configuration=''):
	AS=translate_to_RNA(AS)
	UTR=translate_to_RNA(UTR)
	terminator=translate_to_RNA(terminator)
	CDS=translate_to_RNA(CDS)[:50]
	AS_analysis_dict = dict()
	print features_list
	RNAfoldcall_AS = RNAfold_call(AS+terminator)
	if 'dG_UTR' in features_list:
		AS_analysis_dict['dG_UTR']=dG_UTR
		print 1
	if 'dG_AS' in features_list or 'dG_formAB' in features_list or 'dG_formAA' in features_list:
		dG_AS = extract_dG_ensemble(RNAfoldcall_AS,delete=True)
		AS_analysis_dict['dG_AS']=dG_AS

	if 'dG_AS_AS' in features_list or 'dG_formAA' in features_list:
		dG_AS_AS=extract_dG_ensemble(RNAcofold_call(AS+terminator,AS+terminator),delete=True)
		AS_analysis_dict['dG_AS_AS']=dG_AS_AS
	if 'dG_UTR_AS' in features_list or 'dG_formAB' in features_list or 'rbs_coverage_5' in features_list or 'rbs_coverage_11' in features_list:
		dG_UTR_AS,p_available_in_utr,rbs_coverage_5,rbs_coverage_11 = AS_UTR_analysis(AS+terminator,UTR,p_UTR,CDS,real_RBS)
		AS_analysis_dict['dG_UTR_AS']=dG_UTR_AS
		AS_analysis_dict['p_available_in_utr']=p_available_in_utr
		AS_analysis_dict['rbs_coverage_5']=rbs_coverage_5
		AS_analysis_dict['rbs_coverage_11']=rbs_coverage_11
	if 'dG_formAB' in features_list:
		dG_formAB = dG_UTR_AS - dG_AS - dG_UTR
		AS_analysis_dict['dG_formAB']=dG_formAB
	if 'dG_formAA' in features_list:
		dG_formAA = dG_AS_AS - dG_AS - dG_AS
		AS_analysis_dict['dG_formAA']=dG_formAA
	if 'paired_termini' in features_list:
		paired_termini = calculate_paired_termini(AS,UTR)
		AS_analysis_dict['paired_termini']=paired_termini
	if 'dG_int_seed_avg' in features_list or 'dG_tot_seed_avg' in features_list:
		dG_int_seed_avg,dG_tot_seed_avg,seeds_collection=calculate_seeds(AS+terminator,UTR,CDS,seeds_collection,n=100)
		AS_analysis_dict['dG_int_seed_avg']=dG_int_seed_avg
		AS_analysis_dict['dG_tot_seed_avg']=dG_tot_seed_avg
#		AS_analysis_dict['dG_int_seed_avg']=0
#		AS_analysis_dict['dG_tot_seed_avg']=0

	unique_UTR_id = randomstring(20)
	RNAfold_as_utr = RNAcofold_call(UTR+CDS,AS+terminator,unique_string = unique_UTR_id+'_AS_UTR')
	path_AS_UTR_RNAfold_output = os.getcwd()+'/'+unique_UTR_id+'_AS_UTR'+'_dp.ps'
	p_as_utr = extract_all_nucleotide_binding_probabilities(path_AS_UTR_RNAfold_output,delete=True)

	if UTR+CDS in dict_utr_probs:
		probs_utr = dict_utr_probs[UTR+CDS]
	else:
		RNAfold_utr = RNAfold_call(UTR+CDS,unique_string = unique_UTR_id+'_UTR')
		path_UTR_RNAfold_output = os.getcwd()+'/'+unique_UTR_id+'_UTR'+'_dp.ps'
		probs_utr = extract_all_nucleotide_binding_probabilities(path_UTR_RNAfold_output,delete=True)
		dict_utr_probs[UTR+CDS] = probs_utr
	l_UTR = len(UTR)
	prob_calculations_needed = list(reversed(filter(lambda x:'p_w_as_UTR_' in x, features_list)))+filter(lambda x:'p_w_as_CDS_' in x, features_list)+list(reversed(filter(lambda x:'p_wo_as_UTR_' in x, features_list)))+filter(lambda x:'p_wo_as_CDS_' in x, features_list)
	if len(prob_calculations_needed)>0:
		for pUTR_w in filter(lambda x:'p_w_as_UTR_' in x, prob_calculations_needed):
			nr_UTR = int(pUTR_w.split("UTR_")[1])
			n_nr_UTR = l_UTR+1-nr_UTR
	#		p_utr_nr_UTR = probs_utr[n_nr_UTR]
			p_as_utr_nr_UTR = p_as_utr[n_nr_UTR]
			AS_analysis_dict[pUTR_w] = p_as_utr_nr_UTR
		for pCDS_w in filter(lambda x:'p_w_as_CDS_' in x, prob_calculations_needed):
			nr_CDS = int(pCDS_w.split("CDS_")[1])
			n_nr_CDS = l_UTR+nr_CDS
	#		p_utr_nr_CDS = probs_utr[n_nr_CDS]
			p_as_utr_nr_CDS = p_as_utr[n_nr_CDS]
			AS_analysis_dict[pCDS_w] = p_as_utr_nr_CDS
		for dpUTR_w in filter(lambda x:'dp_w_as_UTR_' in x, prob_calculations_needed):
			nr_UTR = int(dpUTR_w.split("UTR_")[1])
			n_nr_UTR = l_UTR+1-nr_UTR
			p_utr_nr_UTR = probs_utr[n_nr_UTR]
			p_as_utr_nr_UTR = p_as_utr[n_nr_UTR]
			AS_analysis_dict[dpUTR_w] = p_as_utr_nr_UTR-p_utr_nr_UTR
		for dpCDS_w in filter(lambda x:'dp_w_as_CDS_' in x, prob_calculations_needed):
			nr_CDS = int(dpCDS_w.split("CDS_")[1])
			n_nr_CDS = l_UTR+nr_CDS
			p_utr_nr_CDS = probs_utr[n_nr_CDS]
			p_as_utr_nr_CDS = p_as_utr[n_nr_CDS]
			AS_analysis_dict[dpCDS_w] = p_as_utr_nr_CDS-p_utr_nr_CDS
		for pUTR in filter(lambda x:'p_wo_as_UTR_' in x, prob_calculations_needed):
			nr_UTR = int(pUTR.split("UTR_")[1])
			n_nr_UTR = l_UTR+1-nr_UTR
			p_utr_nr_UTR = probs_utr[n_nr_UTR]
	#		p_as_utr_nr_UTR = p_as_utr[n_nr_UTR]
			AS_analysis_dict[pUTR] = p_utr_nr_UTR
		for pCDS in filter(lambda x:'p_wo_as_CDS_' in x, prob_calculations_needed):
			nr_CDS = int(pCDS.split("CDS_")[1])
			n_nr_CDS = l_UTR+nr_CDS
			p_utr_nr_CDS = probs_utr[n_nr_CDS]
	#		p_as_utr_nr_CDS = p_as_utr[n_nr_CDS]
			AS_analysis_dict[pCDS] = p_utr_nr_CDS
#	print AS_analysis_dict
	return AS_analysis_dict,seeds_collection,dict_utr_probs


def RBS_analysis(UTR,rRNA_16S3prime="ACCUCCUUA",delete=True):
	UTR=translate_to_RNA(UTR)
	RBS_part_UTR =  UTR[-18:-4]
	RBS_place_in_UTR_re = re.search(RBS_part_UTR,UTR)
	RBS_place = RBS_place_in_UTR_re.start(0)
	id_RNAcofold_call = randomstring(25)
	RNAcofold_16S_rRNA_RBS_part = RNAcofold_call(RBS_part_UTR,rRNA_16S3prime,unique_string=id_RNAcofold_call)
	path_RNAcofold_output = os.getcwd()+'/'+id_RNAcofold_call+'_dp.ps'
	p_intra_seq1, p_intra_seq2, p_inter_seq1, p_inter_seq2 = extract_all_nucleotide_binding_probabilities(path_RNAcofold_output,RNAcofold = True)
	RBS_center = int(round(sum([p_inter_seq1[x]*x for x in p_inter_seq1.keys()])/sum(p_inter_seq1.values())))
	print UTR+' :  center =   '+str(RBS_center)
	RBS_nucleotides_list = [x+RBS_place-1 for x in range(RBS_center-5+1,RBS_center+5+2)]
	if delete == True:
		try_to_remove_path(path_RNAcofold_output)
	return RBS_nucleotides_list

def calculate_paired_termini(antisense,utr,n=100): #calculates number of paired termini based on n subopt call
	paired_termini_list = []
	RNAsubopt = RNAsubopt_call(antisense,n)
	lAS = len(antisense)
	list_sec_str = []
	for line in RNAsubopt.split('\n'):
	#	print line
		sec_str = re.search(r'([\.\)\(]{'+str(lAS)+'})',line)
		if sec_str:
			list_sec_str.append(sec_str.group(1))
	#calculating paired_termini over 100 samples...
	for subopt_config in list_sec_str: 
		first_config = subopt_config[:lAS/2]
		last_config = subopt_config[lAS/2:]
		open_first = first_config.count('(')
		closed_first = first_config.count(')')
		net_open_first = open_first-closed_first
		open_last = last_config.count('(')
		closed_last = last_config.count(')')
		net_closed_last = closed_last-open_last
		if(net_closed_last == net_open_first):
			paired_termini_list.append(net_closed_last)
	paired_termini = float(sum(paired_termini_list))/float(len(paired_termini_list))
	return(paired_termini)
	
def calculate_seeds(antisense,utr,cds,seeds_collection,n=100): #calculates seeds (int and tot) based on n subopt call
	RNAsubopt = RNAsubopt_call(antisense,n)
	cds = translate_to_RNA(cds)[:50]
	utr = translate_to_RNA(utr)
	lAS = len(antisense)
	list_sec_str = []
	for line in RNAsubopt.split('\n'):
		#print line
		sec_str = re.search(r'([\.\)\(]{'+str(lAS)+'})',line)
		if sec_str:
			list_sec_str.append(sec_str.group(1))
	dG_tot_seed_list = []
	dG_int_seed_list = []
	#seed calculations
	dG_tot_seed, dG_int_seed = 10000,10000
	for subopt_config in list_sec_str: 
	#	print subopt_config
		for l_seed in range(2,6):
			seed_sec_str = re.finditer(r'(?=(\.{'+str(l_seed)+'}))',subopt_config)
			for match in seed_sec_str:
				start_match = match.start(0)
				end_match = match.end(1)
				seed_part = antisense[start_match:end_match]
				if seed_part in seeds_collection.keys():
			#		print seed_part
					dG_int_seed_candidate,dG_tot_seed_candidate = seeds_collection[seed_part]
				else:
					dG_tot_seed_candidate,dG_int_seed_candidate = extract_int_tot_seed(RNAup_call(seed_part,utr+cds))
				if(dG_tot_seed_candidate<dG_tot_seed):
					dG_tot_seed = dG_tot_seed_candidate
				if(dG_int_seed_candidate<dG_int_seed):
					dG_int_seed = dG_int_seed_candidate
					seeds_collection[seed_part] = dG_int_seed_candidate,dG_tot_seed_candidate
	dG_tot_seed_list.append(dG_tot_seed)
	dG_int_seed_list.append(dG_int_seed)
	dG_tot_seed_avg = float(sum(dG_tot_seed_list))/float(len(dG_tot_seed_list))
	dG_int_seed_avg = float(sum(dG_int_seed_list))/float(len(dG_int_seed_list))
	print seeds_collection
	return(dG_int_seed_avg,dG_tot_seed_avg,seeds_collection)
