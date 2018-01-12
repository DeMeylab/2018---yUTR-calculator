import random, time
from functions.general import extract_scale_and_coefficients, coef_scales_from_dict_to_list
from functions.RNA import translate_to_RNA
from functions.riboswitch import analyze_input
from functions.library_tools import recombination_library_lists, degenerate_format_to_list_format,generate_bins,degeneracy_generator,remove_degeneracy,mutate_non_degenerated_part_list, list_to_degenerate_formate, evaluate_degenerated_sequence, evaluate_library_reach_and_resolution
from numpy.random import choice

names_features_utr_analysis = ['dG_EFE','purineAG_in_min3','U_in_min3','A_in_min1','AA_in_min32','CG_in_min32','AC_in_min21','oof_uAUG','GACA_kmer','GG_kmer','CACC_kmer','CA_in_min76','CC_in_min76']
coef_csvfile = '/home/thomas/Dropbox/YeastUTR/coefficients.csv'
scales_csvfile = '/home/thomas/Dropbox/YeastUTR/scales.csv'
coef_dict,scales_dict = extract_scale_and_coefficients(coef_csvfile,scales_csvfile)
coef_list,scales_list,intercept = coef_scales_from_dict_to_list(coef_dict,scales_dict,names_features_utr_analysis)

list_number_of_moves = [1]*8+[2]*4+[3]*2+[4]*1
print(list_number_of_moves)

# Enter CDS
cds = ''

cds = translate_to_RNA(cds)

def library_mutation(pool_candidates,list_available_nucleotides,list_if_nucleotides_are_mandatory,bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept,n_moves=1):
	mutation_moves = ['degeneration']*3+['remove degeneration']*3+['recombination']*6+['point mutations']*6
	candidate_scores = [float(x[0]) for x in pool_candidates.values()]
	tot_prob = sum(candidate_scores)
	prob_candidate_selection_pool = [x/float(tot_prob) for x in candidate_scores] #higher probability for better scoring libraries
	
	candidates_for_mutation = [choice(pool_candidates.keys(),p=prob_candidate_selection_pool),choice(pool_candidates.keys(),p=prob_candidate_selection_pool)]
	for i in range(int(n_moves)):
		move = random.choice(mutation_moves)
		if move == 'degeneration':
			member_mutated = candidates_for_mutation[0]
			member_mutated_list = degenerate_format_to_list_format(member_mutated)
			new_candidate_list=degeneracy_generator(list_available_nucleotides,list_if_nucleotides_are_mandatory,previous_candidate=member_mutated_list)
			new_candidate_degenerated = list_to_degenerate_formate(new_candidate_list)
			while new_candidate_degenerated in pool_candidates.keys():
				
				new_candidate_list=degeneracy_generator(list_available_nucleotides,list_if_nucleotides_are_mandatory,previous_candidate=member_mutated_list)
				new_candidate_degenerated = list_to_degenerate_formate(new_candidate_list)
			inflated_candidate,response_list = evaluate_degenerated_sequence(degenerate_format_to_list_format(new_candidate_degenerated),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
			pool_candidates[new_candidate_degenerated]= evaluate_library_reach_and_resolution(response_list,bins,inflated_candidate)		
			candidates_for_mutation[0] = new_candidate_degenerated
		elif move == 'remove degeneration':
			member_mutated = candidates_for_mutation[0]
			member_mutated_list = degenerate_format_to_list_format(member_mutated)
			if reduce(lambda x, y: x*y, [(len(x)) for x in member_mutated_list]) > 3:
				new_candidate_list=remove_degeneracy(member_mutated_list)
				new_candidate_degenerated = list_to_degenerate_formate(new_candidate_list)
				escape_from_block = 0
				while new_candidate_degenerated in pool_candidates.keys():
					escape_from_block +=  1
					if escape_from_block >5:
						break
					new_candidate_list=remove_degeneracy(member_mutated_list)
					new_candidate_degenerated = list_to_degenerate_formate(new_candidate_list)
				inflated_candidate,response_list = evaluate_degenerated_sequence(degenerate_format_to_list_format(new_candidate_degenerated),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
				pool_candidates[new_candidate_degenerated]= evaluate_library_reach_and_resolution(response_list,bins,inflated_candidate)		
				candidates_for_mutation[0] = new_candidate_degenerated
		elif move == 'point mutations':
			member_mutated = candidates_for_mutation[0]
			member_mutated_list = degenerate_format_to_list_format(member_mutated)
			new_candidate_list=mutate_non_degenerated_part_list(member_mutated_list,list_available_nucleotides)
			new_candidate_degenerated = list_to_degenerate_formate(new_candidate_list)
			escape_from_block = 0
			while new_candidate_degenerated in pool_candidates.keys():
				escape_from_block +=  1
				if escape_from_block >5:
					break
				new_candidate_list=mutate_non_degenerated_part_list(member_mutated_list,list_available_nucleotides)
				new_candidate_degenerated = list_to_degenerate_formate(new_candidate_list)
			inflated_candidate,response_list = evaluate_degenerated_sequence(degenerate_format_to_list_format(new_candidate_degenerated),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
			pool_candidates[new_candidate_degenerated] = evaluate_library_reach_and_resolution(response_list,bins,inflated_candidate)		
			candidates_for_mutation[0] = new_candidate_degenerated
		elif move == 'recombination':
			new_candidate_list1,new_candidate_list2=recombination_library_lists(degenerate_format_to_list_format(candidates_for_mutation[0]),degenerate_format_to_list_format(candidates_for_mutation[1]))
			new_candidate_degenerated1 = list_to_degenerate_formate(new_candidate_list1)
			new_candidate_degenerated2 = list_to_degenerate_formate(new_candidate_list2)
			inflated_candidate1,response_list1 = evaluate_degenerated_sequence(degenerate_format_to_list_format(new_candidate_degenerated1),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
			pool_candidates[new_candidate_degenerated1] = evaluate_library_reach_and_resolution(response_list1,bins,inflated_candidate1)
			inflated_candidate2,response_list2 = evaluate_degenerated_sequence(degenerate_format_to_list_format(new_candidate_degenerated2),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
			pool_candidates[new_candidate_degenerated2] = evaluate_library_reach_and_resolution(response_list2,bins,inflated_candidate2)
			candidates_for_mutation[0] = new_candidate_degenerated1
			candidates_for_mutation[1] = new_candidate_degenerated2
		return pool_candidates	
		
def select_n_best(pool_candidates,n=100):
		pool_candidates_values = [x[0] for x in pool_candidates.values()]
		pool_candidates_values.sort()
		cutoff_value = pool_candidates_values[:n][-1]
		new_pool_candidates = dict()
		for key_candidate in pool_candidates.keys():
			if pool_candidates[key_candidate][0] <= cutoff_value:
				new_pool_candidates[key_candidate] = pool_candidates[key_candidate]
		return new_pool_candidates
#
def coene_library_calculator(list_available_nucleotides,list_if_nucleotides_are_mandatory,bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept,iterations=100,pool_size=100,n_mutations_per_iteration=100,constraint=''):
	pool_candidates = dict()
	for i in range(int(pool_size)):
		new_candidate = degeneracy_generator(list_available_nucleotides,list_if_nucleotides_are_mandatory,first_iteration=True)
		new_candidate=degeneracy_generator(list_available_nucleotides,list_if_nucleotides_are_mandatory,first_iteration=False,previous_candidate=new_candidate)
		new_candidate_degenerated = list_to_degenerate_formate(new_candidate)
		inflated_candidate,response_list = evaluate_degenerated_sequence(degenerate_format_to_list_format(new_candidate_degenerated),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
		pool_candidates[new_candidate_degenerated] = evaluate_library_reach_and_resolution(response_list,bins,inflated_candidate)
	for i in range(int(iterations)):
		print i
		for j in range(int(n_mutations_per_iteration)):
			n_moves_this_mutation = random.choice(list_number_of_moves)
			pool_candidates = library_mutation(pool_candidates,list_available_nucleotides,list_if_nucleotides_are_mandatory,bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept,n_moves=n_moves_this_mutation)
		pool_candidates = select_n_best(pool_candidates,n=pool_size)
	#	print min(pool_candidates.values())
	return pool_candidates
	#	print len(pool_candidates)


#for degenerated_sequence in pool_candidates.keys():
#		print pool_candidates[degenerated_sequence][0]
	#	
	#reduce(lambda x, y: x*y, [(len(x)) for x in previous_candidate])
#	pool_candidates.append([:])
#	print pool_candidates
#print recombination_library_lists(degenerate_format_to_list_format(pool_candidates.keys()[1]),degenerate_format_to_list_format(pool_candidates.keys()[2]))
	
#new_candidate = degeneracy_generator(list_available_nucleotides,list_if_nucleotides_are_mandatory,first_iteration=False,previous_candidate=first_candidate)[:]


#print mutate_non_degenerated_part_list(new_candidate,list_available_nucleotides)

#remove_degeneracy(new_candidate)
utr_list = ['gcatagcaatctaatctaagtttnnnnnnnnnna'] #NO FACULTATIVE NUCLEOTIDES ALLOWED (for the moment)....

minimum = 5 #protein abundance
maximum = 10
n = 8
bins = generate_bins(minimum,maximum,n)
n_iterations = 500
n_pool_size = 100
n_mutations = 100
for utr in utr_list:
	utr = translate_to_RNA(utr)
	list_available_nucleotides,list_if_nucleotides_are_mandatory=analyze_input(utr)
	fh = open('/home/thomas/Dropbox/YeastUTR/coumarate/report_library_creation_'+utr+'_250it_highPA4.txt','w')
	fh.write('n_iterations = '+str(n_iterations)+'\nn_pool_size = '+str(n_pool_size)+'\nn_mutations = '+str(n_mutations)+'\nminimum = '+str(minimum)+'\nmaximum = '+str(maximum)+'\nn_sequences = '+str(n)+'\n')
	pool_candidates = coene_library_calculator(list_available_nucleotides,list_if_nucleotides_are_mandatory,bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept,iterations=n_iterations,pool_size=n_pool_size,n_mutations_per_iteration=n_mutations)
	best_candidates= select_n_best(pool_candidates,n=5)
	candidate_nr = 1
	best_candidates_keys = [x for (y,x) in sorted(zip(best_candidates.values(),best_candidates.keys()))]
	best_candidates_values = [x[0] for x in best_candidates.values()]
	best_candidates_values.sort()
	for candidate in best_candidates_keys:
		score_tuple = best_candidates[candidate]
		fh.write(str(candidate_nr)+') degenerated sequence = '+candidate+' score = '+str(score_tuple[0])+' score_difference = '+str(score_tuple[1])+' score_bins = '+str(score_tuple[2])+' inefficient_bins_coverage_penalty = '+str(score_tuple[3])+'\n')
		inflated_candidate,response_list = evaluate_degenerated_sequence(degenerate_format_to_list_format(candidate),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
		fh.write('sequence,protein_abundance\n')
		for i in range(len(inflated_candidate)):
			fh.write('   -- '+inflated_candidate[i]+','+str(response_list[i])+'\n')
		candidate_nr += 1
 #evaluate_degenerated_sequence(degenerate_format_to_list_format(new_candidate_degenerated),bins,cds,names_features_utr_analysis,coef_list,scales_list,intercept)
