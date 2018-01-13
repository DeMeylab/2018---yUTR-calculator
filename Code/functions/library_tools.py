import random
from RNA import extract_all_nucleotide_binding_probabilities, translate_to_RNA, extract_dG_ensemble, RNAfold_call, re
from general import from_features_to_response_of_PLS
from riboswitch import create_list_with_number_of_replacement


def recombination_library_lists(candidate_list1,candidate_list2):
	degenerated_positions1 = [x for x in range(len(candidate_list1)) if len(candidate_list1[x]) > 1]
	degenerated_positions2 = [x for x in range(len(candidate_list2)) if len(candidate_list2[x]) > 1]
	if len(degenerated_positions1) > 1 and len(degenerated_positions2) > 1 :
		min_recombination = min(degenerated_positions1+degenerated_positions2)
		max_recombination = max(degenerated_positions1+degenerated_positions2)
		place_of_mutation=random.choice(range(min_recombination,max_recombination))
	else:
		place_of_mutation = random.choice(range(1,len(candidate_list1)-1))
	new_candidate_list1 = candidate_list1[:place_of_mutation]+candidate_list2[place_of_mutation:]
	new_candidate_list2 = candidate_list2[:place_of_mutation]+candidate_list1[place_of_mutation:]
	return new_candidate_list1,new_candidate_list2
	
def evaluate_library_reach_and_resolution(response_list,bins_list,inflated_candidate):
	response_list.sort()
#	print response_list
	l_response_list = len(response_list)
	inefficient_bins_coverage_penalty = abs(len(response_list)-float(len(bins_list)-1)) #ideally == 0
#	print response_list
#	print 'eff ='+str(bins_coverage_efficiency)
	if l_response_list == 1:
		objective = 1000
		difference_distance = 1000
		score_bins = 1
	elif l_response_list > 1:
		n_bins_covered = 0
		for i in range(1,len(bins_list)):
			if len([x for x in response_list if bins_list[i-1] < x < bins_list[i]]):
				n_bins_covered +=1
	#	print n_bins_covered
		optimal_distance = float((bins_list[-1]-bins_list[0])/l_response_list)/2
		average_distance = float(sum([response_list[i+1]-response_list[i] for i in range(l_response_list-1)])/(l_response_list-1))
		difference_distance = float(abs(optimal_distance-average_distance)/optimal_distance) # normalized to 1
		n_bins = float(len(bins_list))-float(1)
		not_covered_bins = (n_bins-float(n_bins_covered))/n_bins
#		print 'covered ='+str(covered_bins)
		score_bins = not_covered_bins
		
#		print 'score_bins = '+str(score_bins)
		objective = difference_distance + score_bins*10+inefficient_bins_coverage_penalty
#		print 'n_bins_covered = '+str(n_bins_covered)
#		print 'objective ='+str(objective)
#	print objective,difference_distance,score_bins
	return objective,difference_distance,score_bins,inefficient_bins_coverage_penalty
	
def remove_disfunctional_list(candidate_list):
	new_candidate_list = []
	for i in range(len(candidate_list)):
		if len(candidate_list[i]) == 1 and isinstance(candidate_list[i],list):
			new_candidate_list.append(candidate_list[i][0])
		else:
			new_candidate_list.append(candidate_list[i])
	return new_candidate_list
	
def list_to_degenerate_formate(list_degenerate):
	new_list_degenerate = list_degenerate[:]
	degenerate_dict = {"N":["A","C","G","U"],"V":["A","C","G"],"H":["A","C","U"],"D":["A","G","U"],"B":["C","G","U"],"Y":["C","U"],"R":["A","G"],"K":["G","U"],"M":["A","C"],"S":["C","G"],"W":["A","U"],"A":"A","C":"C","G":"G","U":"U"}
	degenerate_nucleotide_places = [i for i in range(len(list_degenerate)) if len(list_degenerate[i]) > 1]
	for i in degenerate_nucleotide_places:
		for key in degenerate_dict.keys():
			if set(degenerate_dict[key]) == set(list_degenerate[i]):
				new_list_degenerate[i] = key
	new_list_degenerate=remove_disfunctional_list(new_list_degenerate)
	return ''.join(new_list_degenerate)
	
def degenerate_format_to_list_format(degenerate_format):
	new_list_degenerate = []
	degenerate_dict = {"N":["A","C","G","U"],"V":["A","C","G"],"H":["A","C","U"],"D":["A","G","U"],"B":["C","G","U"],"Y":["C","U"],"R":["A","G"],"K":["G","U"],"M":["A","C"],"S":["C","G"],"W":["A","U"],"A":"A","C":"C","G":"G","U":"U"}
	for i in range(len(degenerate_format)):
		new_list_degenerate.append(degenerate_dict[degenerate_format[i]])
	return new_list_degenerate
	
def count_oof_uAUG(utr,start_codon='AUG'):
	counter_oof_uAUG = 0
	interesting_positions = [i for i in [i+1 for i in range(len(utr))] if i % 3 != 0 and i > 2]
	for position in interesting_positions:
		if utr[-position:-(position-3)] == start_codon:
			counter_oof_uAUG+=1
	return(counter_oof_uAUG)
	
def pattern_finder(pattern,utr):
	matched = -1
	match_pattern = '.*'+pattern+'.*'
	if re.match(match_pattern,utr):
		matched = 1
	return(matched)


def generate_bins(minimum,maximum,n):
	bins_list = [float(minimum)]
	for i in range(n):
		added_each_bin = (float(maximum)-float(minimum))/float(n)
		bins_list.append(bins_list[i]+added_each_bin)
	return bins_list

def inflate_list_format_degenerated_sequence(list_formated_degenerated_sequence):
	inflated_list = list_formated_degenerated_sequence[0]
	for i in range(1,len(list_formated_degenerated_sequence)):
		new_inflated_list = []
		for j in range(len(list_formated_degenerated_sequence[i])):
			for k in range(len(inflated_list)):
				new_inflated_list.append(inflated_list[k]+list_formated_degenerated_sequence[i][j])
		inflated_list = new_inflated_list[:]
	return inflated_list	

def remove_degeneracy(previous_candidate):
	n_sequence_list = [x for x in range(len(previous_candidate)) if len(previous_candidate[x])>1]
	if len(n_sequence_list) > 0:
		nucleotide_removing_degeneracy = random.choice(n_sequence_list)
		previous_candidate[nucleotide_removing_degeneracy].remove(random.choice(previous_candidate[nucleotide_removing_degeneracy]))
		if len(previous_candidate[nucleotide_removing_degeneracy]) == 1:
			previous_candidate[nucleotide_removing_degeneracy] = previous_candidate[nucleotide_removing_degeneracy][0]
	else:
		pass
	output_candidate = previous_candidate[:]
	return output_candidate

def mutate_non_degenerated_part_list(previous_candidate,list_available_nucl):
	list_number_replacements = create_list_with_number_of_replacement(list_available_nucl,[1]*len(list_available_nucl))
	output_candidate = previous_candidate[:]
	non_degenerated_nucleotides = [x for x in range(len(previous_candidate)) if len(previous_candidate[x])==1 and len(list_available_nucl[x]) > 1]
	if len(non_degenerated_nucleotides) > 0:
		nucleotide_to_mutate = random.choice(non_degenerated_nucleotides)
		output_candidate[nucleotide_to_mutate] = random.choice([x for x in list_available_nucl[nucleotide_to_mutate] if x != output_candidate[nucleotide_to_mutate]])
	return output_candidate

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
	return [dG_EFE,purineAG_in_min3,U_in_min3,A_in_min1,AA_in_min32,CG_in_min32,AC_in_min21,oof_uAUG,GACA_kmer,GG_kmer,CACC_kmer,CA_in_min76,CC_in_min76]

def degeneracy_generator(list_available_nucleotides,list_if_nucleotides_are_mandatory,first_iteration=False,previous_candidate=''):
	if first_iteration:
		first_candidate_list = []
		for i in range(len(list_available_nucleotides)):
			if list_if_nucleotides_are_mandatory[i] == 1:
				first_candidate_list.append(random.choice(list_available_nucleotides[i]))
			else:
				first_candidate_list.append('')
		return first_candidate_list
	elif previous_candidate != '':
		new_candidate = previous_candidate[:]
		n_sequence = reduce(lambda x, y: x*y, [(len(x)) for x in previous_candidate])
		degenerate_this_nucleotide = random.choice([x for x in range(len(previous_candidate)) if set(previous_candidate[x]) != set(list_available_nucleotides[x]) ])
		possible_replacements = [x for x in list_available_nucleotides[degenerate_this_nucleotide] if x not in previous_candidate[degenerate_this_nucleotide]]
		if len(previous_candidate[degenerate_this_nucleotide])==1:
				
			new_candidate[degenerate_this_nucleotide] = [previous_candidate[degenerate_this_nucleotide][0],random.choice(possible_replacements)]
		else:
			new_candidate[degenerate_this_nucleotide] = new_candidate[degenerate_this_nucleotide]+[random.choice(possible_replacements)]
		return new_candidate
	
def evaluate_degenerated_sequence(candidate,bins_list,cds,names_features_utr_analysis,coef_list,scales_list,intercept):
	inflated_candidate = inflate_list_format_degenerated_sequence(candidate)
	response_list = []
	for candidate_utr_seq in inflated_candidate:
		utr_features = analyze_UTR(candidate_utr_seq,cds)
		response = from_features_to_response_of_PLS(utr_features,names_features_utr_analysis,coef_list,scales_list,intercept)
		response_list.append(response)
	return inflated_candidate,response_list
