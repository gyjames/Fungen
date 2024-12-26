#!/data/software/conda/bin/python


import networkx as nx
import multiprocessing as mp
import os, sys
import random
from collections import deque, Counter
import parasail
import re
from time import time


def sort_best_seq(read_id_lst, read_seq_lst):
	seq_lst = [(read_id_lst[read_i], read_seq_lst[read_i]) for read_i in range(len(read_id_lst))]
	seq_lst.sort(key = lambda x:len(x[1]), reverse = True)
	seq_lst_len = len(seq_lst)
	kmer_dic = {}
	for seq_name, seq in seq_lst[: min(10 ** 3, len(seq_lst))]:
		for k in range(len(seq) - 10 + 1):
			kmer = seq[k : k + 10]
			if kmer not in kmer_dic:
				kmer_dic[kmer] = 0
			kmer_dic[kmer] += 1
	seq_lst_scored = []
	seq_number = 0
	for seq_name, seq in seq_lst:
		seq_number += 1
		if seq_number < min(10 ** 3, len(seq_lst)):
			seq_score = sum([kmer_dic[seq[k : k + 10]] for k in range(len(seq) - 10 + 1) if seq[k : k + 10] in kmer_dic])
			seq_lst_scored.append((seq_name, seq, seq_score))
		else:
			seq_score = len(seq)
			seq_lst_scored.append((seq_name, seq, seq_score))
	seq_lst_scored.sort(key = lambda x:x[2], reverse = True)
	best_seq_name, best_seq, _ = seq_lst_scored[0]
	read_id_lst_sorted = [read_id for read_id, read_seq, seq_score in seq_lst_scored]
	read_seq_lst_sorted = [read_seq for read_id, read_seq, seq_score in seq_lst_scored]
	return read_id_lst_sorted, read_seq_lst_sorted

def two_seq_similarity(sequence1, sequence2):
	aaset = ''.join([i for i in list(set(sequence1+sequence2))])
	aamatrix = parasail.matrix_create(aaset,2,-2)
	parasail_result = parasail.sg_trace_scan_16(sequence1,sequence2,5,1,aamatrix)
	parasail_cigar = str(parasail_result.cigar.decode,'utf-8')
	result = re.split(r'[=DXSMI]+', parasail_cigar)
	cigar_result,i = [],0
	for cigar_num in result[:-1]:
		i+=len(cigar_num)+1
		cigar_result.append((int(cigar_num),parasail_cigar[i-1]))
	#print(cigar_result)
	nt_match,nt_mismatch = 0,0
	while True:
		if cigar_result[0][-1] != '=':
			cigar_result.remove(cigar_result[0])
		if len(cigar_result) == 0:
			break
		if cigar_result[-1][-1] != '=':
			cigar_result.remove(cigar_result[-1])
		if len(cigar_result) == 0:
			break
		if cigar_result[0][-1] == '=' and cigar_result[-1][-1] == '=':
			break
	#print(cigar_result)
	if len(cigar_result) == 0:
		ll = 0
	else:
		cigar_result.sort(key=lambda x:x[0], reverse = True)
		ll = 0
		for l,m in cigar_result:
			if m == '=':
				ll = l
				break
	sequence2_in_align = 0
	cig_i = 0
	for cig in cigar_result:
		cig_i += 1
		if cig[1] == '=':
			nt_match += cig[0]
		else:
			if cig[0] < 10:
				nt_mismatch += cig[0]
						

		if cig[1] == '=' or cig[1] == 'I':
			sequence2_in_align += cig[0]
	if nt_match+nt_mismatch == 0:
		return 0
	#if nt_match < len(sequence1)/3 or nt_match < len(sequence2)/3:
		#return 0
	if nt_match < 20:
		return 0
	
	#print(nt_match/(nt_match+nt_mismatch), nt_match, nt_match+nt_mismatch)
	return nt_match/(nt_match+nt_mismatch)
def two_seq_similarity_equal_max(sequence1, sequence2):
	aaset = ''.join([i for i in list(set(sequence1+sequence2))])
	aamatrix = parasail.matrix_create(aaset,2,-2)
	parasail_result = parasail.sg_trace_scan_16(sequence1,sequence2,5,1,aamatrix)
	parasail_cigar = str(parasail_result.cigar.decode,'utf-8')
	result = re.split(r'[=DXSMI]+', parasail_cigar)
	cigar_result,i = [],0
	for cigar_num in result[:-1]:
		i+=len(cigar_num)+1
		cigar_result.append((int(cigar_num),parasail_cigar[i-1]))
	#print(cigar_result)
	nt_match,nt_mismatch = 0,0
	while True:
		if cigar_result[0][-1] != '=':
			cigar_result.remove(cigar_result[0])
		if len(cigar_result) == 0:
			break
		if cigar_result[-1][-1] != '=':
			cigar_result.remove(cigar_result[-1])
		if len(cigar_result) == 0:
			break
		if cigar_result[0][-1] == '=' and cigar_result[-1][-1] == '=':
			break
	#print(cigar_result)
	if len(cigar_result) == 0:
		ll = 0
	else:
		cigar_result.sort(key=lambda x:x[0], reverse = True)
		ll = 0
		for l,m in cigar_result:
			if m == '=':
				ll = l
				break
	sequence2_in_align = 0
	cig_i = 0
	longest_qual = 0
	for cig in cigar_result:
		cig_i += 1
		if cig[1] == '=':
			if cig[0] > longest_qual:
				longest_qual = cig[0]
			nt_match += cig[0]
		else:
			if cig[0] < 10:
				nt_mismatch += cig[0]


		if cig[1] == '=' or cig[1] == 'I':
			sequence2_in_align += cig[0]
	if nt_match+nt_mismatch == 0:
		return 0, 0
	#if nt_match < len(sequence1)/3 or nt_match < len(sequence2)/3:
		#return 0
	if nt_match < 20:
		return 0, 0

	#print(nt_match/(nt_match+nt_mismatch), nt_match, nt_match+nt_mismatch)
	return nt_match/(nt_match+nt_mismatch), longest_qual

def score_sequence(seq_lst):
	seq_lst.sort(key = lambda x:len(x[1]), reverse = True)
	kmer_dic = {}
	for seq_name, seq in seq_lst[: min(10 ** 4, len(seq_lst))]:
		for k in range(len(seq) - 10 + 1):
			kmer = seq[k : k + 10]
			if kmer not in kmer_dic:
				kmer_dic[kmer] = 0
			kmer_dic[kmer] += 1
	seq_lst_scored = []
	seq_number = 0
	for seq_name, seq in seq_lst:
		seq_number += 1
		if seq_number < min(10 ** 4, len(seq_lst)):
			seq_score = sum([kmer_dic[seq[k : k + 10]] for k in range(len(seq) - 10 + 1)])
			seq_lst_scored.append((seq_name, seq, seq_score))
		else:
			seq_lst_scored.append((seq_name, seq, len(seq_lst) - seq_name))
	seq_lst_scored.sort(key = lambda x:x[2], reverse = True)
	return seq_lst_scored


def add_seq_to_matrix(ref_seq, query_seq):
	aaset = ''.join([i for i in list(set(query_seq + ref_seq))])
	aamatrix = parasail.matrix_create(aaset, 2, -2)
	parasail_result = parasail.sg_trace_scan_16(query_seq, ref_seq, 5, 1, aamatrix)
	parasail_cigar = str(parasail_result.cigar.decode, 'utf-8')
	#print(parasail_cigar)
	result = re.split(r'[=DXSMI]+', parasail_cigar)
	i, cigar_tuples = 0, []
	for length in result[:-1]:
		i += len(length)
		type_ = parasail_cigar[i]
		i += 1
		cigar_tuples.append((int(length), type_ ))
	r_index, q_index, q_aln, r_aln,  = 0, 0, [], []
	query_aln_lst = []
	cigar_tuple_len = len(cigar_tuples)
	query_site_i = 0
	start_I = ''

	match_num = 0
	aln_num = 0
	for cigar_n in range(len(cigar_tuples)):
		length_, type_ = cigar_tuples[cigar_n]
		if type_ == 'D':
			if cigar_n !=0 and cigar_n != cigar_tuple_len - 1:
				aln_num += length_
				for k in range(length_):
					query_aln_lst.append('-')
			else:
				for k in range(length_):
					query_aln_lst.append('|')
		if type_ == '=' or type_ == 'X':
			aln_num += length_
			if type_ == '=':
				match_num += length_
			for k in range(length_):
				query_aln_lst.append(query_seq[query_site_i])
				query_site_i += 1
		if type_ == 'I':
			if cigar_n !=0:
				if cigar_n != cigar_tuple_len - 1:
					aln_num += length_
				for k in range(length_):
					if query_aln_lst[-1] != '-' and query_aln_lst[-1] != '|':
						query_aln_lst[-1] += query_seq[query_site_i]
						query_site_i += 1
					else:
						query_aln_lst[-1] = query_seq[query_site_i]
						query_site_i += 1
			if cigar_n == 0:
				for k in range(length_):
					start_I += query_seq[query_site_i]
					query_site_i += 1
	if start_I != '':
		if query_aln_lst[0] != '-' and query_aln_lst[0] != '|':
			query_aln_lst[0] = start_I + query_aln_lst[0]
		else:
			query_aln_lst[0] = start_I
	#print((match_num + 1) / (aln_num + 1))
	if (match_num + 1) / (aln_num + 1) < 0.7:
		return False, ['|' for i in range(len(ref_seq))]
	if (match_num + 1) / (aln_num + 1) >= 0.7:
		return True, query_aln_lst

def identify_sv_sites(seq_aln_matrix):
	read_lst_len = len(seq_aln_matrix)
	sv_lst = []
	seq_len = len(seq_aln_matrix[0])
	#print([len(seq_lst) for seq_lst in seq_aln_matrix])
	for site_lo in range(seq_len):
		sv_site_lst = []
		site_lst = [seq_lst[site_lo] for seq_lst in seq_aln_matrix]
		site_lst_clean = [w for w in site_lst if w != '|']
		site_count_lst = list(Counter(site_lst_clean).most_common())
		if len(site_count_lst) <= 1 or len(site_lst_clean) < len(site_lst)/10:
			continue
		site_count_lst = [(site_nt, site_num) for site_nt, site_num in site_count_lst if site_num > 1]
		if len(site_count_lst) < 1:
			continue
		if len(site_count_lst) == 1:
			continue
		for site_nt, site_num in site_count_lst:
			same_kmer_lst = []
			sv_nt_lst = []
			for read_num in range(read_lst_len):
				if site_lst[read_num] == site_nt:
					sv_nt_lst.append(read_num)
					sv_kmer = ''.join(seq_aln_matrix[read_num][max(0, site_lo - 3) : min(site_lo + 4, seq_len)])
					if '|' not in sv_kmer:
						same_kmer_lst.append(sv_kmer)
			if len(same_kmer_lst) > 1:
				same_kmer_count = Counter(same_kmer_lst).most_common(1)
				if same_kmer_count[0][1] > max(1, len(site_lst_clean) * 0.05):#more than 2 kmers sopport and abundance > 5 %
					sv_site_lst.append((site_nt, sv_nt_lst))
		if len(sv_site_lst) > 1:
			sv_lst.append((site_lo, sv_site_lst))
		#if len(sv_lst) > 1:
		#	print(site_lo, 'sv_lst readnum', [(site, len(read_lst)) for site, read_lst in sv_lst], len(site_lst))
	return sv_lst

def extract_heaviest_path(sv_lst):
	connected_reads = []
	sv_lst_len = len(sv_lst)
	for sv_site_i in range(sv_lst_len - 1):
		site_lo1, site_sv_lst1 = sv_lst[sv_site_i]
		site_lo2, site_sv_lst2 = sv_lst[sv_site_i + 1]
		connected_reads_new = []
		for nt_site2, read_lst2 in site_sv_lst2:
			site2_path = []
			for nt_site1, read_lst1 in site_sv_lst1:
				if site_lo1 == '|' or site_lo2 == '|':
					continue
				read_lst_ttl = read_lst1 + read_lst2
				overlap_num = len(read_lst_ttl) - len(list(set(read_lst_ttl)))
				if sv_site_i == 0:
					site2_path.append([[(site_lo1, nt_site1), (site_lo2, nt_site2)], overlap_num])
				else:
					for cr_num in range(len(connected_reads)):
						if connected_reads[cr_num][0][-1] == (site_lo1, nt_site1):
							path_extend = list(connected_reads[cr_num][0])
							path_extend.append((site_lo2, nt_site2))
							path_weight = connected_reads[cr_num][1] + overlap_num
							site2_path.append([path_extend, path_weight])
			if site2_path != []:
				site2_path.sort(key = lambda x:x[-1], reverse = True)
				connected_reads_new.append(site2_path[0])
		connected_reads = connected_reads_new
	connected_reads.sort(key = lambda x:x[-1], reverse = True)
	if connected_reads != []:
		heaviest_path = connected_reads[0][0]
	else:
		heaviest_path = []
	return heaviest_path

def purify_corelated_sites(seq_aln_matrix, corelated_sites):
	purified_corelated_sites = []
	t_read_num = len(seq_aln_matrix)
	aln_len = len(seq_aln_matrix[0])
	for corelated_sites_part in corelated_sites:
		corelated_sites_part_lst = list(corelated_sites_part)
		corelated_sites_part_lst.sort(key = lambda x:x[0])
		#print(corelated_sites_part_lst)
		main_species_count = Counter()
		for site_i, site_nt, read_id_lst in corelated_sites_part:
			main_species_count.update(read_id_lst)
		core_reads = [read_id for read_id, read_support_num in list(main_species_count.most_common()) if read_support_num > max(2, len(corelated_sites_part) * 0.7)]
		#print('core_reads', core_reads)
		site_support_core_reads = []
		for site_i, site_nt, read_id_lst in corelated_sites_part_lst:
			read_lst_core_site = core_reads + list(read_id_lst)
			overlapped_read_num = len(read_lst_core_site) - len(list(set(read_lst_core_site)))
			#print(overlapped_read_num, len(core_reads), len(read_id_lst), overlapped_read_num > len(core_reads) * 0.3, overlapped_read_num > len(read_id_lst) * 0.5)
			if overlapped_read_num > max(2, len(core_reads) * 0.7) and overlapped_read_num > max(2, len(read_id_lst) * 0.5):
				if site_support_core_reads != []:
					if site_i > site_support_core_reads[-1][0] + 6:
						site_support_core_reads.append((site_i, site_nt, overlapped_read_num))
					elif site_i == site_support_core_reads[-1][0]:
						if overlapped_read_num > site_support_core_reads[-1][-1]:
							site_support_core_reads[-1] = (site_i, site_nt, overlapped_read_num)
				elif site_support_core_reads == []:
					site_support_core_reads.append((site_i, site_nt, overlapped_read_num))
		effecitive_site_num = len([site_nt for site_i, site_nt, overlapped_read_num in site_support_core_reads if site_nt != '-'])
		#print('effecitive_site_num', aln_len * 0.001, effecitive_site_num, site_support_core_reads)
		if effecitive_site_num > max(3, aln_len * 0.001):
			site_support_core_reads.sort(key = lambda x:x[0])
			purified_corelated_sites.append([(site_i, site_nt) for site_i, site_nt, overlap_read_num in site_support_core_reads if site_nt != '-'])
	#print('purified_corelated_sites', purified_corelated_sites)
	return purified_corelated_sites


def extract_corelated_sites(seq_aln_matrix, sv_lst, heaviest_path):
	sv_num = len(sv_lst)
	#print('sv_num', sv_num)
	sv_num_half = int(sv_num/2)
	#sv_lst = sv_lst[max(0, sv_num_half - 100): min(sv_num_half + 100, sv_num)]#use as much as 200 sv sites
	heaviest_path_set = set(heaviest_path)
	corelated_sites = []
	for sv_site_i in range(len(sv_lst)):
		for sv_site_j in range(len(sv_lst)):
			if sv_site_j <= sv_site_i:
				continue
			site_lo1, site_sv_lst1 = sv_lst[sv_site_i]
			site_lo2, site_sv_lst2 = sv_lst[sv_site_j]
			for nt_site1, read_lst1 in site_sv_lst1:
				for nt_site2, read_lst2 in site_sv_lst2:
					read_site1 = (site_lo1, nt_site1)
					read_site2 = (site_lo2, nt_site2)
					if read_site1 in heaviest_path_set or read_site2 in heaviest_path_set:
						continue
					if site_lo2 - site_lo1 <= 6:
						continue
					overlapped_read_num = len(read_lst1 + read_lst2) - len(list(set(read_lst1 + read_lst2)))
					if overlapped_read_num >= int(0.3 * max(len(read_lst1), len(read_lst2))) + 1 and overlapped_read_num >= int(0.9 * min(len(read_lst1), len(read_lst2))) + 1:
						best_match_read_lst2, best_match_num2 = [], 0 #make sure lst1 and lst2 are the best matches in both sites
						for nt_site, read_lst in site_sv_lst2:
							match_num = len(read_lst + read_lst1) - len(set(list(read_lst + read_lst1)))
							if match_num > best_match_num2:
								best_match_read_lst2, best_match_num2 = read_lst, match_num

						best_match_read_lst1, best_match_num1 = [], 0
						for nt_site, read_lst in site_sv_lst1:
							match_num = len(read_lst + read_lst2) - len(list(set(read_lst + read_lst2)))
							if match_num > best_match_num1:
								best_match_read_lst1, best_match_num1 = read_lst, match_num
						if best_match_read_lst2 == read_lst2 and best_match_read_lst1 == read_lst1:
							corelated_sites.append(((site_lo1, nt_site1, tuple(read_lst1)), (site_lo2, nt_site2, tuple(read_lst2))))
	#print(corelated_sites) #Get corelated sites and read_lst
	
	combine_bool = True
	while combine_bool == True:
		combine_bool = False
		corelated_sites_new = []
		for read_site_set in corelated_sites:
			hit_bool = False
			for n in range(len(corelated_sites_new)):
				hit_bool = False
				for read_site in list(read_site_set):
					if read_site in corelated_sites_new[n]:
						corelated_sites_new[n].update(read_site_set)
						hit_bool = True
						break
				if hit_bool:
					combine_bool = True
					break
			if not hit_bool:
				corelated_sites_new.append(set(read_site_set))
		corelated_sites = corelated_sites_new
	#print(corelated_sites)
	purified_corelated_sites_lst = purify_corelated_sites(seq_aln_matrix, corelated_sites)
	#for corelated_sites in purified_corelated_sites_lst:
	#	print(corelated_sites)
	return purified_corelated_sites_lst

def add_site_to_path(heaviest_path, corelated_sites):
	corelated_sites_dic = {site_i: site_nt for site_i, site_nt in corelated_sites}
	corelated_heaviest_path = []
	for site_i, site_nt in heaviest_path:
		if site_i in corelated_sites_dic:
			corelated_heaviest_path.append((site_i, corelated_sites_dic[site_i]))
		elif site_i not in corelated_sites_dic:
			corelated_heaviest_path.append((site_i, site_nt))
	return corelated_heaviest_path


def extract_max_flow(seq_aln_matrix, sv_lst):
	heaviest_path_lst = [] 
	heaviest_path = extract_heaviest_path(sv_lst)
	heaviest_path_lst.append(heaviest_path)
	#print(heaviest_path)
	corelated_sites = extract_corelated_sites(seq_aln_matrix, sv_lst, heaviest_path)
	for site_lst in corelated_sites:
		heaviest_path_lst.append(add_site_to_path(heaviest_path, site_lst))
	return heaviest_path_lst


def get_consensus_seq(seq_aln_matrix, sv_lst):
	sv_consensus_dic = {}
	for site_i, site_nt_lst in sv_lst:
		site_nt_lst.sort(key = lambda x:len(x[1]), reverse = True)
		sv_consensus_dic[site_i] = site_nt_lst[0][0]
	seq_consensus = []
	for site_lo in range(len(seq_aln_matrix[0])):
		if site_lo in sv_consensus_dic:
			site_nt_most_freq = sv_consensus_dic[site_lo]
			seq_consensus.append(site_nt_most_freq)
		if site_lo not in sv_consensus_dic:
			site_nt_lst = [seq_aln[site_lo] for seq_aln in seq_aln_matrix if seq_aln[site_lo] != '|']
			if len(site_nt_lst) == 0:
				site_nt_most_freq = '-'
			else:
				site_nt_most_freq = list(Counter(site_nt_lst).most_common())[0][0]
			seq_consensus.append(site_nt_most_freq)
	return seq_consensus


def aln_seq_to_consensus(seq_aln_lst, path_consensus):
	aln_score = 0
	tt_num = 0
	for site_lo in range(len(seq_aln_lst)):
		if seq_aln_lst[site_lo] == path_consensus[site_lo]:
			aln_score += 1
		if seq_aln_lst[site_lo] != '|':
			tt_num += 1
	if aln_score < 5:
		return 0
	return (aln_score + 1) / (tt_num + 1)


def matrix_seq(seq_id_lst, seq_lst):
	aln_seq_id_lst = [seq_id_lst[0]]
	ref_seq = seq_lst[0]
	#print('ref seq len', len(ref_seq))
	seq_len = len(ref_seq)
	seq_aln_matrix = [[w for w in ref_seq]]
	start_time = time()
	for query_seq_i in range(1, len(seq_lst)):
		query_seq = seq_lst[query_seq_i]
		query_id = seq_id_lst[query_seq_i]
		add_bool, query_seq_aln_lst = add_seq_to_matrix(ref_seq, query_seq)
		if add_bool:
			#print(len(query_seq_aln_lst))
			aln_seq_id_lst.append(query_id)
			seq_aln_matrix.append(query_seq_aln_lst)
	#print(len(seq_aln_matrix), len(aln_seq_id_lst))
	return aln_seq_id_lst, seq_aln_matrix


def classify_reads(heaviest_path_lst, aln_seq_id_lst, seq_aln_matrix):
	path_num = len(heaviest_path_lst)
	path_site_num = len(heaviest_path_lst[0])
	informative_site_lst = [[] for i in range(path_num)]
	for path_site_i in range(path_site_num):
		if len(list(set([heaviest_path_lst[path_i][path_site_i] for path_i in range(path_num)]))) != 1:
			for path_i in range(path_num):
				informative_site_lst[path_i].append(heaviest_path_lst[path_i][path_site_i])
	#print(informative_site_lst)
	info_len = len(informative_site_lst[0])

	path_reads_dic = {i:[] for i in range(path_num)}
	for read_i in range(len(seq_aln_matrix)):
		#print('read_i', read_i)
		read_id = aln_seq_id_lst[read_i]
		read_aln_seq_lst = seq_aln_matrix[read_i]
		read_path_score = []
		for path_i in range(path_num):
			#print('path_i', path_i)
			effc_num = 0
			score_i = 0
			for path_site_i in range(info_len):
				site_i, site_nt = informative_site_lst[path_i][path_site_i]
				site_in_read = read_aln_seq_lst[site_i]
				if site_in_read != '|':
					effc_num += 1
				if site_in_read == site_nt:
					score_i += 1
			#print(score_i, effc_num)
			#if (score_i + 1) / (effc_num + 1) > 0.8 and score_i >= 2:
			read_path_score.append((path_i, score_i, effc_num))
		read_path_score.sort(key = lambda x:x[1], reverse = True)
		#print('read_path_score', read_path_score)
		if read_path_score != []:
			path_reads_dic[read_path_score[0][0]].append(read_i)
	cluster_include_ref = []
	cluster_include_ref_aln = []
	cluster_include_ref_ord = []
	cluster_heaviest_path = []
	#print('path_reads_dic', path_reads_dic)
	path_reads_lst = [path_reads_dic[path_i] for  _, path_i in enumerate(path_reads_dic)]
	path_reads_lst.sort(key = lambda x:len(x), reverse = True)
	path_reads_ord = path_reads_lst[0]
	cluster_include_ref_ord = path_reads_ord
	cluster_include_ref = [aln_seq_id_lst[read_i] for read_i in path_reads_dic[path_i]]
	cluster_include_ref_aln = [seq_aln_matrix[read_i] for read_i in path_reads_dic[path_i]]
	cluster_heaviest_path = heaviest_path_lst[path_i]
			
	return cluster_include_ref, cluster_include_ref_aln, cluster_include_ref_ord, cluster_heaviest_path

def polish_sequencing_error(cluster_include_ref, cluster_include_ref_aln):
	sv_lst = identify_sv_sites(cluster_include_ref_aln[: min(200, len(cluster_include_ref))])
	seq_consensus = get_consensus_seq(cluster_include_ref_aln[: min(200, len(cluster_include_ref))], sv_lst)
	heaviest_path = extract_heaviest_path(sv_lst)
	heaviest_path_dic = {site_i: site_nt for site_i, site_nt in heaviest_path}
	corrected_seq_consensus = []
	for site_i in range(len(seq_consensus)):#Correct seq_consensus using heaviest_path
		site_nt = seq_consensus[site_i]
		if site_i in heaviest_path_dic:
			corrected_seq_consensus.append(heaviest_path_dic[site_i])
			#print(site_nt, heaviest_path_dic[site_i])
		else:
			corrected_seq_consensus.append(site_nt)
	corrected_reads_include_ref = []
	for read_i in range(len(cluster_include_ref)):
		read_id_i = cluster_include_ref[read_i]
		read_seq_aln = cluster_include_ref_aln[read_i]
		corrected_read_aln = []
		aln_num = 0
		match_num = 0
		for site_i in range(len(read_seq_aln)):
			if read_seq_aln[site_i] != '|':
				aln_num += 1
				if corrected_seq_consensus[site_i] == read_seq_aln[site_i]:
					match_num += 1
				if corrected_seq_consensus[site_i] != '-':
					corrected_read_aln.append(corrected_seq_consensus[site_i])
		#print(match_num, aln_num)
		if (match_num + 1) / (aln_num + 1) > 0.7:
			corrected_reads_include_ref.append((read_id_i, ''.join(corrected_read_aln)))
	#print(corrected_reads_include_ref)
	return corrected_reads_include_ref

def recombine_corrected_read_cluster(corrected_read_cluster_lst):
	corrected_read_cluster_lst_recombine = []
	for read_cluster in corrected_read_cluster_lst:
		#print(read_cluster[0][1])
		best_hit_cluster, best_similarity, longest_equal = 0, 0, 0
		for recombined_read_cluster_i in range(len(corrected_read_cluster_lst_recombine)):
			recombined_read_cluster = corrected_read_cluster_lst_recombine[recombined_read_cluster_i]
			similarity, equal_max = two_seq_similarity_equal_max(read_cluster[0][1], recombined_read_cluster[0][1])
			if similarity > best_similarity:
				best_similarity = similarity
				best_hit_cluster = recombined_read_cluster_i
				longest_equal = equal_max
		if best_similarity > 0.99 or longest_equal > 100:
			corrected_read_cluster_lst_recombine[best_hit_cluster] += read_cluster
		else:
			corrected_read_cluster_lst_recombine.append(read_cluster)
		#print(best_similarity)
	return corrected_read_cluster_lst_recombine


def compare_represents_mp_bak(rep_info):
	cluster_n, repseq_n, cluster_m, repseq_m = rep_info
	similarity, equal_max = two_seq_similarity_equal_max(repseq_n, repseq_m)
	return cluster_n, cluster_m, similarity, equal_max

def compare_represents_mp(rep_info):
	cluster_n, cluster_m = rep_info
	repseq_n = CORRECTED_READ_CLUSTER_LST_P[cluster_n]
	repseq_m = CORRECTED_READ_CLUSTER_LST_P[cluster_m]
	similarity, equal_max = two_seq_similarity_equal_max(repseq_n, repseq_m)
	return cluster_n, cluster_m, similarity, equal_max

def initializer_recombine_corrected_read_cluster_mp(corrected_read_cluster_lst_p):
	global CORRECTED_READ_CLUSTER_LST_P
	CORRECTED_READ_CLUSTER_LST_P = corrected_read_cluster_lst_p

def recombine_corrected_read_cluster_mp(corrected_read_cluster_lst, threads):
	corrected_read_cluster_lst_recombine = []
	compare_info_mp = []
	corrected_read_cluster_lst_p = [r[0][1] for r in corrected_read_cluster_lst]
	compare_dic = {}
	for n in range(len(corrected_read_cluster_lst)):
		for m in range(len(corrected_read_cluster_lst)):
			#compare_info_mp.append((n, corrected_read_cluster_lst[n][0][1], m, corrected_read_cluster_lst[m][0][1]))
			compare_info_mp.append((n, m))
			if len(compare_info_mp) == 10 ** 7 or (n + 1 == len(corrected_read_cluster_lst) and m + 1 == len(corrected_read_cluster_lst)):
				p = mp.Pool(processes = threads, initializer = initializer_recombine_corrected_read_cluster_mp(corrected_read_cluster_lst_p))
				rep_compare_results = p.map(compare_represents_mp, compare_info_mp)
				p.close()
				p.join()

				for cluster_n, cluster_m, similarity, equal_max in rep_compare_results:
					if similarity <= 0.99 and equal_max <= 100:
						continue
					if cluster_n not in compare_dic:
						compare_dic[cluster_n] = {}
					compare_dic[cluster_n][cluster_m] = (similarity, equal_max)
				compare_info_mp = []
	raw_i_dic = {}
	for read_cluster_i in range(len(corrected_read_cluster_lst)):
		read_cluster = corrected_read_cluster_lst[read_cluster_i]
		#print(read_cluster[0][1])
		best_hit_cluster, best_similarity, longest_equal = 0, 0, 0
		for recombined_read_cluster_i in range(len(corrected_read_cluster_lst_recombine)):
			raw_recombined_read_cluster_i = raw_i_dic[recombined_read_cluster_i]
			if read_cluster_i not in compare_dic:
				continue
			if raw_recombined_read_cluster_i not in compare_dic[read_cluster_i]:
				continue
			similarity, equal_max = compare_dic[read_cluster_i][raw_recombined_read_cluster_i]
			if similarity > best_similarity:
				best_similarity = similarity
				best_hit_cluster = recombined_read_cluster_i
				longest_equal = equal_max
		if best_similarity > 0.99 or longest_equal > 100:
			corrected_read_cluster_lst_recombine[best_hit_cluster] += read_cluster
		else:
			corrected_read_cluster_lst_recombine.append(read_cluster)
			raw_i_dic[len(corrected_read_cluster_lst_recombine) - 1] = read_cluster_i
		#print(best_similarity)
	return corrected_read_cluster_lst_recombine
						

def correct_seqs_singleMP(seq_info):
	start_time = time()
	seq_id_lst, seq_lst = seq_info
	corrected_read_cluster_lst = []
	correct_bool = False
	while len(seq_id_lst) >= 3:#Max flow to extract and correct the reads belong to heaviest_path
		aln_seq_id_lst, seq_aln_matrix = matrix_seq(seq_id_lst, seq_lst)
		sv_lst = identify_sv_sites(seq_aln_matrix)
		#print(sv_lst)
		heaviest_path_lst = extract_max_flow(seq_aln_matrix, sv_lst)
		corrected_read_clust = classify_reads(heaviest_path_lst, aln_seq_id_lst, seq_aln_matrix)
		#print('corrected_read_clust', corrected_read_clust)
		seq_consensus = get_consensus_seq(seq_aln_matrix, sv_lst)
		cluster_include_ref, cluster_include_ref_aln, cluster_include_ref_ord, cluster_heaviest_path = classify_reads(heaviest_path_lst, aln_seq_id_lst, seq_aln_matrix)
		if len(cluster_include_ref) >= 3:
			polished_read_lst = polish_sequencing_error(cluster_include_ref, cluster_include_ref_aln)
			corrected_read_cluster_lst.append(polished_read_lst)
		elif len(cluster_include_ref) < 3:
			polished_read_lst = []
		#print(cluster_include_ref, cluster_include_ref_ord)
		if len(polished_read_lst) <= 3:
			break
		polished_read_set = set([read_id for read_id,seq in polished_read_lst])
		seq_lst = [seq_lst[read_i] for read_i in range(len(seq_id_lst)) if seq_id_lst[read_i] not in polished_read_set]
		seq_id_lst = [seq_id_lst[read_i] for read_i in range(len(seq_id_lst)) if seq_id_lst[read_i] not in polished_read_set]
		if len(seq_id_lst) > 0:
			seq_id_lst, seq_lst = sort_best_seq(seq_id_lst, seq_lst)
	corrected_read_cluster_lst_recombine = recombine_corrected_read_cluster(corrected_read_cluster_lst)
	#print(time() - start_time)
	return corrected_read_cluster_lst_recombine



def correct_seqs_singleT(seq_info):
	seq_id_lst, seq_lst = seq_info
	corrected_read_cluster_lst = []
	while len(seq_id_lst) >= 3:#Max flow to extract and correct the reads belong to heaviest_path
		#seq_chunk_num = int(len(seq_id_lst) / threads) + 1
		seq_chunk_num = 10000
		max_threads = int(((len(seq_id_lst) - 1)) / seq_chunk_num)
		if len(seq_id_lst) - max_threads * seq_chunk_num < 5000:
			max_threads -= 1
		max_threads = max(max_threads, 1)
		for thread_i in range(max_threads):
			if thread_i != max_threads - 1:
				seq_id_lst_chunk = seq_id_lst[thread_i * seq_chunk_num: min((thread_i + 1) * seq_chunk_num, len(seq_id_lst))]
				seq_lst_chunk = seq_lst[thread_i * seq_chunk_num: min((thread_i + 1) * seq_chunk_num, len(seq_lst))]
				se_id_lst_chunk, seq_lst_chunk = sort_best_seq(seq_id_lst_chunk, seq_lst_chunk)
			elif thread_i == max_threads - 1:
				seq_id_lst_chunk = seq_id_lst[thread_i * seq_chunk_num: len(seq_id_lst)]
				seq_lst_chunk = seq_lst[thread_i * seq_chunk_num: len(seq_lst)]
				seq_id_lst_chunk, seq_lst_chunk = sort_best_seq(seq_id_lst_chunk, seq_lst_chunk)
			seq_info = (seq_id_lst_chunk, seq_lst_chunk)
			seq_correct_results_lst = correct_seqs_singleMP(seq_info)
			correct_read_id_set = set([])
			corrected_bool = False
			#print(seq_correct_results_lst)
			for corrected_read_cluster in seq_correct_results_lst:
				if len(corrected_read_cluster) >=3:
					corrected_bool = True
					corrected_read_cluster_lst.append(corrected_read_cluster)
					correct_read_id_set.update([read_id for read_id, seq in corrected_read_cluster])
		seq_lst = [seq_lst[read_i] for read_i in range(len(seq_id_lst)) if seq_id_lst[read_i] not in correct_read_id_set]
		seq_id_lst = [read_id for read_id in seq_id_lst if read_id not in correct_read_id_set]
	#	print('#', len(corrected_read_cluster_lst), corrected_bool)
		if not corrected_bool:
			break
	#	print('remain', len(seq_id_lst))
	corrected_read_cluster_lst_recombine = recombine_corrected_read_cluster(corrected_read_cluster_lst)
	return corrected_read_cluster_lst_recombine



def correct_seqs_MP(seq_id_lst, seq_lst, threads):
	corrected_read_cluster_lst = []
	while len(seq_id_lst) >= 3:#Max flow to extract and correct the reads belong to heaviest_path
		#seq_chunk_num = int(len(seq_id_lst) / threads) + 1
		#print(len(seq_id_lst))
		seq_chunk_num = 10000
		max_threads = int(((len(seq_id_lst) - 1)) / seq_chunk_num)
		if len(seq_id_lst) - max_threads * seq_chunk_num < 5000:
			max_threads -= 1
		max_threads = max(max_threads, 1)
		seq_mp_info = []
		for thread_i in range(max_threads):
			if thread_i != max_threads - 1:
				seq_id_lst_chunk = seq_id_lst[thread_i * seq_chunk_num: min((thread_i + 1) * seq_chunk_num, len(seq_id_lst))]
				seq_lst_chunk = seq_lst[thread_i * seq_chunk_num: min((thread_i + 1) * seq_chunk_num, len(seq_lst))]
				seq_id_lst_chunk, seq_lst_chunk = sort_best_seq(seq_id_lst_chunk, seq_lst_chunk)
				seq_mp_info.append((seq_id_lst_chunk, seq_lst_chunk))
			elif thread_i == max_threads - 1:
				seq_id_lst_chunk = seq_id_lst[thread_i * seq_chunk_num: len(seq_id_lst)]
				seq_lst_chunk = seq_lst[thread_i * seq_chunk_num: len(seq_lst)]
				seq_id_lst_chunk, seq_lst_chunk = sort_best_seq(seq_id_lst_chunk, seq_lst_chunk)
				seq_mp_info.append((seq_id_lst_chunk, seq_lst_chunk))
		p = mp.Pool(processes = min(threads, max_threads))
		seq_correct_results_lst = p.map(correct_seqs_singleMP, seq_mp_info)
		p.close()
		p.join()
		correct_read_id_set = set([])
		corrected_bool = False
		for corrected_read_cluster_lst_chunk in seq_correct_results_lst:
			for corrected_read_cluster in corrected_read_cluster_lst_chunk:
				if len(corrected_read_cluster) >=3:
					corrected_bool = True
					corrected_read_cluster_lst.append(corrected_read_cluster)
					correct_read_id_set.update([read_id for read_id, seq in corrected_read_cluster])
		if not corrected_bool:
			break
		seq_lst = [seq_lst[read_i] for read_i in range(len(seq_id_lst)) if seq_id_lst[read_i] not in correct_read_id_set]
		seq_id_lst = [read_id for read_id in seq_id_lst if read_id not in correct_read_id_set]
		#print('remain', len(seq_id_lst))
	#print(len(corrected_read_cluster_lst))
	#print('combine')
	combine_time = time()
	#corrected_read_cluster_lst_recombine = recombine_corrected_read_cluster(corrected_read_cluster_lst)
	corrected_read_cluster_lst_recombine = recombine_corrected_read_cluster_mp(corrected_read_cluster_lst, threads)
	#print('combine_finish', time() - combine_time)
	return corrected_read_cluster_lst_recombine


def get_read_lst(outdir):
	read_lst_lst = []
	with open('{}/ref_seqs_3st.txt'.format(outdir)) as read_clusters:
		for line in read_clusters:
			read_lst = []
			read_id_lst = line.split()
			for read_id in read_id_lst:
				read_lst.append(int(read_id))
			read_lst_lst.append(read_lst)
	return read_lst_lst
				

def extract_reads_from_files(read_lst_lst,  chunk_read_num):
	reads_clust_lst = []
	reads_file_dic = {}
	for lst_i in range(len(read_lst_lst)):
		read_lst = read_lst_lst[lst_i]
		reads_clust_lst.append(['' for i in range(len(read_lst))])
		for read_i in range(len(read_lst)):
			seq_num_raw = read_lst[read_i]
			seq_num = abs(seq_num_raw)
			file_name, record_num = int((seq_num -1) / chunk_read_num), (seq_num - 1) % chunk_read_num + 1
			if file_name not in reads_file_dic:
				reads_file_dic[file_name] = []
			reads_file_dic[file_name].append((record_num, seq_num_raw, lst_i, read_i))


	tt = 0
	reads_file_lst = []
	for i,file_name in enumerate(reads_file_dic):
		reads_file_lst.append((file_name, reads_file_dic[file_name]))
	random.shuffle(reads_file_lst)

	for file_name,read_record_lst in reads_file_lst:
		#print('scanned read file {}'.format(file_name))
		read_record_lst.sort(key = lambda x:x[0])
		tt = 0
		read_record_deque = deque(read_record_lst)
		if len(read_record_deque) != 0:
			fst_read_num, fst_seq_num_raw, fst_lst_i, fst_seq_i = read_record_deque.popleft()
		else:
			break
		line_num = 0
		with open('{}/seq_record/{}'.format(OUTDIR, file_name)) as seq_record_file:
			for seq_record_line in seq_record_file:
				line_num +=1
				if line_num == fst_read_num:
					seq, seq_rev = seq_record_line.split('\t')[2], seq_record_line.split('\t')[4].strip()
					if fst_seq_num_raw > 0:
						reads_clust_lst[fst_lst_i][fst_seq_i] = seq
					if fst_seq_num_raw < 0:
						reads_clust_lst[fst_lst_i][fst_seq_i] = seq_rev
					if len(read_record_deque) == 0:
						break
					else:
						fst_read_num, fst_seq_num_raw, fst_lst_i, fst_seq_i = read_record_deque.popleft()

	#seq_chunk_data.sort(key = lambda x:x[1])
	return reads_clust_lst



def correct_reads(args, outdir):
	start_time = time()
	print("\033[31m**************\033[0m", 'Correct reads', "\033[31m**************\033[0m")
	global OUTDIR
	OUTDIR = outdir
	threads = args.t
	read_chunk_size = args.c
	read_lst_lst = get_read_lst(outdir)
	#read_lst_lst = read_lst_lst[0 : 3]#use for test
	corrected_read_cluster_lst = []
	read_lst_lst_chunk = []
	read_num_chunk = 0
	chunk_num = 0
	for cluster_n in range(len(read_lst_lst)):
		read_num_chunk += len(read_lst_lst[cluster_n])
		read_lst_lst_chunk.append(read_lst_lst[cluster_n])
		if read_num_chunk >= 10 ** 6 or cluster_n + 1 == len(read_lst_lst) or len(read_lst_lst_chunk) >= 10 ** 5:
			#print('current_chunk', cluster_n)
			reads_file_lst_chunk = extract_reads_from_files(read_lst_lst_chunk,  read_chunk_size)
			chunk_num += 1
			read_info_mp = []
			#print(len(reads_file_lst), reads_file_lst[-3], read_lst_lst[-3])
			for cluster_i in range(len(reads_file_lst_chunk)):
				read_id_lst, seq_lst = read_lst_lst_chunk[cluster_i], reads_file_lst_chunk[cluster_i]
				complex_score = len(read_id_lst) * len(seq_lst[0])
				#print(complex_score)
				if complex_score < 100000 * 100:
				#if True:
					read_info_mp.append((complex_score, read_id_lst, seq_lst))
				else:
					cluster_start_time = time()
					for corrected_read_lst in correct_seqs_MP(read_id_lst, seq_lst, threads):
						corrected_read_cluster_lst.append(corrected_read_lst)
					print("\033[31m[Correct large clusters]\033[0m", 'corrected cluster {}, time used: {}  s'.format(cluster_i + cluster_n, time() - cluster_start_time))
			del read_lst_lst_chunk, reads_file_lst_chunk
			read_info_mp.sort(key = lambda x:x[0], reverse = True)
			read_info_mp = [(read_id_lst, seq_lst) for complex_score, read_id_lst, seq_lst in read_info_mp]
			print("\033[31m[Correct small clusters]\033[0m", 'correcting chunk {}'.format(chunk_num))
			cluster_start_time = time()
			p = mp.Pool(processes = threads)
			correct_clusters_mp_results = p.map(correct_seqs_singleT, read_info_mp)
			p.close()
			p.join()
			for correct_clusters in correct_clusters_mp_results:
				for corrected_read_lst in correct_clusters :
					corrected_read_cluster_lst.append(corrected_read_lst)
			del read_info_mp, correct_clusters_mp_results
			print("\033[31m[Correct small clusters]\033[0m", 'finished correct small clusters, time used: {}  s'.format(time() - cluster_start_time))
			read_num_chunk = 0
			read_lst_lst_chunk = []
			reads_file_lst_chunk = []	
	corrected_read_cluster_lst.sort(key = lambda x:len(x), reverse = True)
	read_name_index = {}
	with open('{}/name_num_record.txt'.format(outdir)) as read_name_record_data:
		for line in read_name_record_data:
			line_spt = line.split()
			read_id, read_name = int(line_spt[0]), line_spt[1]
			read_name_index[read_id] = read_name
			
	cluster_results = []
	correct_results = []
	cluster_represent_seqs = []
	for cluster_i in range(len(corrected_read_cluster_lst)):
		for read_id, read_seq in corrected_read_cluster_lst[cluster_i]:
			if read_id > 0:
				read_name = read_name_index[read_id]
			else:
				read_name = ''.join((read_name_index[-read_id], '|Corr_Reverse'))
			cluster_results.append('{}\t{}\n'.format(cluster_i, read_name))
			correct_results.append('@{}\n{}\n+\n{}\n'.format(read_name, read_seq, 'F' * len(read_seq)))
			if len(cluster_represent_seqs) <= cluster_i:
				cluster_represent_seqs.append('>{}_Cluster-{}\n{}\n'.format(read_name, cluster_i, read_seq))
	cluster_results_file = open('{}/cluster_results.txt'.format(outdir), 'w')
	cluster_results_file.writelines(cluster_results)
	cluster_results_file.close()

	correct_results_file = open('{}/corrected_seqs.fastq'.format(outdir), 'w')
	correct_results_file.writelines(correct_results)
	correct_results_file.close()

	cluster_represent_seqs_file = open('{}/cluster_represents.fasta'.format(outdir), 'w')
	cluster_represent_seqs_file.writelines(cluster_represent_seqs)
	cluster_represent_seqs_file.close()
	print("\033[31m[Correct finished]\033[0m", 'corrected {} clusters, time used: {}  s'.format(len(cluster_represent_seqs), time() - start_time))
