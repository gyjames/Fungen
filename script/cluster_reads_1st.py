#!/data/software/conda/bin/python
import multiprocessing as mp
import os,sys
from time import time
from script.get_minimizers import compress_sequence, get_minimizers, rev_com_seq, rev_string
import math
from collections import Counter
import gzip
import parasail
import re


def generate_minimizers(seq_info):
	seq_num, seq_name_raw, seq, k, w = seq_info
	compressed_seq = compress_sequence(seq)
	compressed_seq_rev = rev_com_seq(compressed_seq)
	query_minimizers_raw = [minimizer for minimizer, i in get_minimizers(compressed_seq, k, w)]
	query_minimizers = [int(minimizer) for minimizer in query_minimizers_raw]
	query_minimizers_rev = [int(minimizer) for minimizer, i in get_minimizers(compressed_seq_rev,k, w)]
	return(seq_num, seq_name_raw, query_minimizers, query_minimizers_rev)

def generate_record_seq(seq_name_record_info):
	seq_name_record, _ = seq_name_record_info
	return seq_name_record, _




def longest_mini_mer(query_minimizers, ref_minimizers):
	ref_minimizer_dic = {}
	ref_minimizer_i = 0
	for ref_minimizer in ref_minimizers:
		try: ref_minimizer_dic[ref_minimizer].append(ref_minimizer_i)
		except:  ref_minimizer_dic[ref_minimizer] = [ref_minimizer_i]
		ref_minimizer_i += 1
	same_region = []
	last_query_end = 0
	for query_minimizer_i in range(len(query_minimizers)):
		#print('query_minimizer_i, last_query_end', query_minimizer_i, last_query_end)
		query_minimizer = query_minimizers[query_minimizer_i]
		if query_minimizer_i < last_query_end:
			continue
		recent_same_region = ((0, 0),(0,0))
		if query_minimizer in ref_minimizer_dic:
			#print(ref_minimizer_dic[query_minimizer])
			#recent_same_region = ((0, 0),(0,0))
			#(ref_region_start, ref_region_end), (query_region_start, query_region_end) = recent_same_region
			for ref_region_start in ref_minimizer_dic[query_minimizer]:
				query_region_start, query_region_end = query_minimizer_i, query_minimizer_i
				ref_region_end = ref_region_start
				while query_region_start - 1 >= 0 and ref_region_start - 1 >= 0 and query_minimizers[query_region_start - 1] == ref_minimizers[ref_region_start - 1]:
					ref_region_start -= 1
					query_region_start -= 1
				while ref_region_end + 1 < len(ref_minimizers) and query_region_end + 1 < len(query_minimizers) and query_minimizers[query_region_end + 1] == ref_minimizers[ref_region_end + 1]:
					ref_region_end += 1
					query_region_end += 1
				#print((query_region_start, query_region_end), (ref_region_start, ref_region_end))
				if query_region_end - query_region_start == ref_region_end - ref_region_start:
					if query_region_end - query_region_start > recent_same_region[0][1] - recent_same_region[0][0] and query_region_end - query_region_start >= 3:
						recent_same_region = ((query_region_start, query_region_end), (ref_region_start, ref_region_end))
		if recent_same_region != ((0, 0),(0,0)):
			if recent_same_region not in same_region:
				same_region.append(recent_same_region)
				last_query_end = recent_same_region[0][1]
	minimizer_match = []
	#print('same_region', same_region)
	for ((query_region_start, query_region_end), (ref_region_start, ref_region_end)) in same_region:
		for i in range(query_region_end - query_region_start + 1):
			#minimizer_match.append((query_minimizers[query_region_start + i], query_region_start + i), (ref_minimizers[ref_region_start + i], ref_region_start + i))
			minimizer_match.append((query_region_start + i, ref_region_start + i))
	return minimizer_match
			


def generate_2mer(minimizer1, minimizer2):
	mini_2mer = int(minimizer1) * (10**6) + int(minimizer2)
	#print(minimizer1, minimizer2, mini_2mer)
	return mini_2mer

def generate_3mer(minimizer1, minimizer2, minimizer3):
	#mini_3mer = int(minimizer1) * (10**12) + int(minimizer2) * (10**6) + int(minimizer3)
	#mini_3mer = ''.join((minimizer1, '.', minimizer2, '.', minimizer3))
	#mini_3mer = int(minimizer1) + int(minimizer2) + int(minimizer3)
	mini_3mer = int(''.join((minimizer1, '0' * max(0, 6 - len(minimizer2)),  minimizer2, '0' * max(0, 6 - len(minimizer3)), minimizer3)))
	#print(minimizer1, minimizer2, mini_2mer)
	return mini_3mer

def decord_minimizers(record_data):
	num_data, minimizer_data = record_data.split('\t')[:2]
	read_num = int(num_data)
	read_numv = -int(num_data)
	query_minimizer_data = minimizer_data.split('|')[0].strip()
	query_minimizer_data_rev = minimizer_data.split('|')[1].strip()
	#query_minimizers = [int(w) for w in query_minimizer_data.split()]
	#query_minimizers_rev = [int(w) for w in query_minimizer_data_rev.split()]
	query_minimizers = [w for w in query_minimizer_data.split()]
	query_minimizers_rev = [w for w in query_minimizer_data_rev.split()]
	return read_num, read_numv, query_minimizers, query_minimizers_rev
	



def compare_file2file(minimizer_record_file1, minimizer_record_file2, outdir):
	compare_start_time = time()
	minimizer_list1 = open('{}/minimizer_record_sort/{}'.format(outdir, minimizer_record_file1)).readlines()
	minimizer_list2 = open('{}/minimizer_record_sort/{}'.format(outdir, minimizer_record_file2)).readlines()
	#seq_list1 = open('{}/seq_record_sort/{}'.format(minimizer_record_file1)).readlines()
	#seq_list2 = open('{}/seq_record_sort/{}'.format(minimizer_record_file2)).readlines()
	
	seq_dic1 = {}
	read_num_1 = 0
	with open('{}/seq_record_sort/{}'.format(outdir, minimizer_record_file1)) as seq_record_data:
		for seq_line in seq_record_data:
			read_num_1 += 1
			seq, seq_v = seq_line.split()[1], seq_line.split()[2].strip()
			seq_dic1[read_num_1] = seq
			seq_dic1[-read_num_1] = seq_v

	seq_dic2 = {}
	read_num_2 = 0
	with open('{}/seq_record_sort/{}'.format(outdir, minimizer_record_file2)) as seq_record_data:
		for seq_line in seq_record_data:
			read_num_2 += 1
			seq, seq_v = seq_line.split()[1], seq_line.split()[2].strip()
			seq_dic2[read_num_2] = seq
			seq_dic2[-read_num_2] = seq_v

	mini_mer_index = {}
	ref_minimizer_list_list = [[]]
	ref_name_list = ['']
	#print('minimizer_lst_finish')
	read_num_1 = 1
	mini_mer_num_dic2 = {}
	for data2_line in minimizer_list2:
		ref_name, ref_name_v, ref_minimizers_lst, ref_minimizers_lstv = decord_minimizers(data2_line)
		for n in range(len(ref_minimizers_lst) - 2):
			minimizer1, minimizer2, minimizer3 = ref_minimizers_lst[n], ref_minimizers_lst[n + 1], ref_minimizers_lst[n + 2]
			ref_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
			if ref_mini_3mer not in mini_mer_num_dic2:
				mini_mer_num_dic2[ref_mini_3mer] = 0
			mini_mer_num_dic2[ref_mini_3mer] += 1
	for data1_line in minimizer_list1:
		ref_name, ref_name_v, ref_minimizers_lst, ref_minimizers_lstv = decord_minimizers(data1_line)
		ref_minimizer_list_list.append(ref_minimizers_lst)
		ref_name_list.append(ref_name)
		for n in range(len(ref_minimizers_lst) - 2):
			minimizer1, minimizer2, minimizer3 = ref_minimizers_lst[n], ref_minimizers_lst[n + 1], ref_minimizers_lst[n + 2]
			ref_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
			if ref_mini_3mer not in mini_mer_index:
				mini_mer_index[ref_mini_3mer] = []
			mini_mer_index[ref_mini_3mer].append(read_num_1)
		read_num_1 += 1
	discard_mini_mer = [(0, 0) for i in range(5000)]
	dele_ref_minimizers = []
	
	tt_3mer_num  = 1
	for i, ref_mini_3mer in enumerate(mini_mer_index):
		mer_num2 = 0
		if ref_mini_3mer in mini_mer_num_dic2:
			mer_num2 = mini_mer_num_dic2[ref_mini_3mer]
		mer_num = len(mini_mer_index[ref_mini_3mer]) + mer_num2
		if mer_num == 1:
			dele_ref_minimizers.append(ref_mini_3mer)
			continue
		tt_3mer_num += 1
		if mer_num > discard_mini_mer[-1][1]:
			discard_mini_mer[-1] = (ref_mini_3mer, mer_num)
			discard_mini_mer.sort(key = lambda x:x[1], reverse = True)

	discard_mini_mer = discard_mini_mer[: min([max(min(200, int(tt_3mer_num / (10 ** 3))), 20), int(tt_3mer_num / (10))])]
	#print('discard_mini_mer', discard_mini_mer)
	discard_mini_mer_set = set([ref_mini_3mer for ref_mini_3mer,mer_num in discard_mini_mer if mer_num > len(minimizer_list1) / 100])
	for ref_mini_3mer in dele_ref_minimizers:
		del mini_mer_index[ref_mini_3mer]
	for ref_mini_3mer in list(discard_mini_mer_set):
		del mini_mer_index[ref_mini_3mer]
	mini_mer_num_dic2 = list([])
	del mini_mer_num_dic2
	
	line_num = 0
	match_output = []
	output_clust_data = []
	for line_num2 in range(len(minimizer_list2)):
		data2_line_num = line_num2 + 1
		data2_line = minimizer_list2[line_num2]
		keep_bool = True
		match_result = []
		output_data_line = []
		line_num += 1
		match_time = time()
		query_name, query_namev, query_minimizers_lst, query_minimizers_lstv = decord_minimizers(data2_line)

		match_count = Counter()
		query_minimizers_lst_len = len(query_minimizers_lst)
		mini_3mer_lst = []
		for m in range(query_minimizers_lst_len - 2):
			minimizer1, minimizer2, minimizer3 = query_minimizers_lst[m], query_minimizers_lst[m + 1], query_minimizers_lst[m + 2]
			query_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
			mini_3mer_lst.append(query_mini_3mer)

		discard_mer_num = 0
		for query_mini_3mer in list(set(mini_3mer_lst)):
			if query_mini_3mer in discard_mini_mer_set:
				discard_mer_num += 1
				continue
			if query_mini_3mer in mini_mer_index:
				match_count.update(list(set(mini_mer_index[query_mini_3mer])))
		match_count_lst = list(match_count.most_common())
		match_num_cutoff = 1 + int(0.1 * (len(list(set(mini_3mer_lst))) - discard_mer_num))

		query_seq = seq_dic2[line_num]
		query_7mer_set = set([query_seq[i : i + 7] for i in range(len(query_seq) - 7 + 1)])

		for match_result in match_count_lst:
			#print(match_result)
			ref_record_num, ref_match_num = match_result
			if ref_match_num < match_num_cutoff:
				break
			ref_seq = seq_dic1[ref_record_num]
			dif_transcript_ratio = confirm_not_different_transcript(ref_seq, query_7mer_set)
			if dif_transcript_ratio > 0.75:
				#print('dif_transcript_ratio, ref_match_num', dif_transcript_ratio, ref_match_num, match_num_cutoff, query_name, ref_name_list[ref_record_num])
				output_data_line.append('{}|{}|{}\t'.format(ref_name_list[ref_record_num], query_name, ref_match_num))
			#print('{}|{}|{}\t'.format(ref_name_list[ref_record_num], query_name, ref_match_num))

		match_count = Counter()
		query_minimizers_lstv_len = len(query_minimizers_lstv)
		counter_time = time()
		mini_3mer_lst = []
		for m in range(query_minimizers_lstv_len - 2):
			minimizer1, minimizer2, minimizer3 = query_minimizers_lstv[m], query_minimizers_lstv[m + 1], query_minimizers_lstv[m + 2]
			query_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
			mini_3mer_lst.append(query_mini_3mer)
		discard_mer_num = 0
		for query_mini_3mer in list(set(mini_3mer_lst)):
			if query_mini_3mer in discard_mini_mer_set:
				discard_mer_num += 1
				continue
			if query_mini_3mer in mini_mer_index:
				match_count.update(list(set(mini_mer_index[query_mini_3mer])))
		match_count_lstv = list(match_count.most_common())
		match_num_cutoffv = 1 + int(0.1 * (len(list(set(mini_3mer_lst))) - discard_mer_num))

		query_seqv = seq_dic2[-line_num]
		query_7mer_setv = set([query_seq[i : i + 7] for i in range(len(query_seqv) - 7 + 1)])

		for match_resultv in match_count_lstv:
			ref_record_numv, ref_match_numv = match_resultv
			if ref_match_numv < match_num_cutoffv:
				break
			ref_seqv = seq_dic1[ref_record_numv]
			dif_transcript_ratiov = confirm_not_different_transcript(ref_seqv, query_7mer_setv)
			if dif_transcript_ratiov > 0.75:
				#print('dif_transcript_ratiov, ref_match_num', dif_transcript_ratiov, ref_match_numv, match_num_cutoffv, len(list(set(mini_3mer_lst))), query_name, ref_name_list[ref_record_numv])
				output_data_line.append('{}|{}|{}\t'.format(-ref_name_list[ref_record_numv], query_name, ref_match_numv))
		if output_data_line != []:
			output_clust_data.append(''.join((''.join(output_data_line), '\n')))

	clust_record_file = open('{}/match_record/{}-{}'.format(outdir, minimizer_record_file1, minimizer_record_file2), 'w')
	clust_record_file.writelines(output_clust_data)
	clust_record_file.close()
	print("\033[31m[Cluster sequence]\033[0m", 'finished compare {} to {}'.format(minimizer_record_file1, minimizer_record_file2), 'time_used: %.3f s'%( time() - compare_start_time))



def find_best_reads(reads_in_cluster_raw):
	#start_time = time()
	read_set = set([])
	reads_in_cluster = []
	for read_info in reads_in_cluster_raw:
		read_line_num = abs(read_info[0])
		if read_line_num not in read_set:
			reads_in_cluster.append(read_info)
			read_set.add(read_line_num)
	#print('refine s1', time() - start_time)
	if len(reads_in_cluster) == 1:
		return reads_in_cluster[0][0], []
	if len(reads_in_cluster) == 2:
		#read_name_in_cluster1, read_line1, seq_line1, mer_num1 = reads_in_cluster[0]
		#read_name_in_cluster2, read_line2, seq_line2, mer_num2 = reads_in_cluster[1]
		read_name_in_cluster1, read_line1, mer_num1 = reads_in_cluster[0]
		read_name_in_cluster2, read_line2, mer_num2 = reads_in_cluster[1]
		best_read = read_name_in_cluster1
		read_name_in_cluster, mer_num = read_name_in_cluster2, mer_num2
		_, _, ref_minimizers_lst1, _ = decord_minimizers(read_line1)
		_, _, ref_minimizers_lst2, _ = decord_minimizers(read_line2)
		if len(ref_minimizers_lst2) > len(ref_minimizers_lst1):
			best_read = read_name_in_cluster2
			read_name_in_cluster, mer_num = read_name_in_cluster1, mer_num1
		reads_in_cluster_refine = [(read_name_in_cluster, mer_num)]
		return best_read, reads_in_cluster_refine

	read_chunk = []
	min_minimizer_lst_len = 10 ** 4
	#for read_name_in_cluster, read_line, seq_line, mer_num in reads_in_cluster[: min(100, len(reads_in_cluster))]:
	for read_name_in_cluster, read_line,  mer_num in reads_in_cluster[: max(min(100, len(reads_in_cluster)), 20)]:
		ref_name, ref_name_v, ref_minimizers_lst, ref_minimizers_lstv = decord_minimizers(read_line) 
		if len(ref_minimizers_lst) < min_minimizer_lst_len:
			min_minimizer_lst_len = len(ref_minimizers_lst)
		if read_name_in_cluster > 0:
			ref_3mer_lst = []
			for n in range(len(ref_minimizers_lst) - 2):
				minimizer1, minimizer2, minimizer3 = ref_minimizers_lst[n], ref_minimizers_lst[n + 1], ref_minimizers_lst[n + 2]
				ref_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
				ref_3mer_lst.append(ref_mini_3mer)
			#read_chunk.append((read_name_in_cluster, ref_3mer_lst, seq_line.split()[0], mer_num))
			read_chunk.append((read_name_in_cluster, ref_3mer_lst, mer_num))
		if read_name_in_cluster < 0:
			ref_3mer_lst = []
			for n in range(len(ref_minimizers_lstv) - 2):
				minimizer1, minimizer2, minimizer3 = ref_minimizers_lstv[n], ref_minimizers_lstv[n + 1], ref_minimizers_lstv[n + 2]
				ref_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
				ref_3mer_lst.append(ref_mini_3mer)
			#read_chunk.append((read_name_in_cluster, ref_3mer_lst, seq_line.split()[1], mer_num))
			read_chunk.append((read_name_in_cluster, ref_3mer_lst,  mer_num))
	mini_3mer_num_dic = {}
	#for read_name_in_cluster, ref_3mer_lst, read_seq, mer_num in read_chunk:
	for read_name_in_cluster, ref_3mer_lst, mer_num in read_chunk:
		for mini_3mer in list(set(ref_3mer_lst)):
			if mini_3mer in mini_3mer_num_dic:
				mini_3mer_num_dic[mini_3mer] += 1
			else:
				mini_3mer_num_dic[mini_3mer] = 1
	mini_3mer_num_lst = [(mini_3mer, mini_3mer_num_dic[mini_3mer]) for _, mini_3mer in enumerate(mini_3mer_num_dic)]
	mini_3mer_num_lst.sort(key = lambda x:x[1], reverse = True)
	min_trust_mini_3mer_num =  min(int(min_minimizer_lst_len / 5), len(mini_3mer_num_lst))
	mini_3mer_trust_lst = []
	mini_3mer_cutoff = max(min(100, len(reads_in_cluster)) / 10, 2)
	for mini_3mer, mini_3mer_num in mini_3mer_num_lst:
		if mini_3mer_num < mini_3mer_cutoff and len(mini_3mer_trust_lst) > min_trust_mini_3mer_num:
			break
		mini_3mer_trust_lst.append(mini_3mer)
	trust_mini_3mer_set = set(mini_3mer_trust_lst)	
	#print('len(trust_mini_3mer_set), min_minimizer_lst_len', len(list(trust_mini_3mer_set)), min_minimizer_lst_len, len(mini_3mer_num_lst))
	read_chunk_match_lst = []
	#for read_name_in_cluster, ref_3mer_lst, ref_seq, mer_num in read_chunk:
	for read_name_in_cluster, ref_3mer_lst, mer_num in read_chunk:
		read_match_score = 0
		for mini_3mer in ref_3mer_lst:
			if mini_3mer in trust_mini_3mer_set:
				read_match_score += 1
		#read_chunk_match_lst.append((read_name_in_cluster, read_match_score, ref_seq))
		read_chunk_match_lst.append((read_name_in_cluster, read_match_score))
	read_chunk_match_lst.sort(key = lambda x:x[1], reverse = True)
	#print(read_chunk_match_lst)
	#best_read, best_read_seq = read_chunk_match_lst[0][0], read_chunk_match_lst[0][2]
	best_read = read_chunk_match_lst[0][0]
	reads_in_cluster_refine = []
	#for read_name_in_cluster, read_line, seq_line, mer_num in reads_in_cluster:
	for read_name_in_cluster, read_line, mer_num in reads_in_cluster:
		if abs(read_name_in_cluster) == abs(best_read):
			continue
		reads_in_cluster_refine.append((read_name_in_cluster, mer_num))
	#print('refine s2', time() - start_time)
	return best_read, reads_in_cluster_refine
	
			
def confirm_not_different_transcript(ref_seq, query_7mer_set):
	ref_7mer_lst = [ref_seq[i : i + 7] for i in range(len(ref_seq) - 7 + 1)]
	ref_7mer_in_num = 0
	for kmer in ref_7mer_lst:
		if kmer in query_7mer_set:
			ref_7mer_in_num += 1
	return ref_7mer_in_num / (len(ref_7mer_lst) + 1)
	

def compare_file_self(minimizer_record_file1, r_n, outdir):
	#print(minimizer_record_file1)
	compare_start_time = time()
	minimizer_list = open('{}/minimizer_record_sort/{}'.format(outdir, minimizer_record_file1)).readlines()
	seq_dic = {}
	read_num_1 = 0
	seq_lst = open('{}/seq_record_sort/{}'.format(outdir, minimizer_record_file1)).readlines()
	for seq_line in seq_lst:
		read_num_1 += 1
		seq, seq_v = seq_line.split()[1], seq_line.split()[2].strip()
		seq_dic[read_num_1] = seq
		seq_dic[-read_num_1] = seq_v
	mini_mer_index = {}
	ref_name_list = ['']
	#print('minimizer_lst_finish')
	read_num_2 = 0
	clust_results_dic = {}
	clustered_reads = set([])
	#Make mini-3mer index
	for data_line in minimizer_list:
		read_num_2 += 1
		ref_name, ref_name_v, ref_minimizers_lst, ref_minimizers_lstv = decord_minimizers(data_line)
		ref_name_list.append(ref_name)
		for n in range(len(ref_minimizers_lst) - 2):
			minimizer1, minimizer2, minimizer3 = ref_minimizers_lst[n], ref_minimizers_lst[n + 1], ref_minimizers_lst[n + 2]
			ref_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
			if ref_mini_3mer not in mini_mer_index:
				mini_mer_index[ref_mini_3mer] = []
			mini_mer_index[ref_mini_3mer].append(read_num_2)
	#tt_3mer_num = len(list(mini_mer_index))
	tt_3mer_num = 1
	discard_mini_mer = [(0, 0) for i in range(500)]
	dele_ref_minimizers = []
	for i, ref_mini_3mer in enumerate(mini_mer_index):
		mer_num = len(mini_mer_index[ref_mini_3mer])
		if mer_num == 1:
			dele_ref_minimizers.append(ref_mini_3mer)
			continue
		tt_3mer_num += 1
		if mer_num > discard_mini_mer[-1][1]:
			discard_mini_mer[-1] = (ref_mini_3mer, mer_num)
			discard_mini_mer.sort(key = lambda x:x[1], reverse = True)
	for ref_mini_3mer in dele_ref_minimizers:
		del mini_mer_index[ref_mini_3mer]
	#print('index 3mer finished', tt_3mer_num)
	discard_mini_mer = discard_mini_mer[: min([max(min(200, int(tt_3mer_num / (10 ** 3))), 20), int(tt_3mer_num / (10))])]
	#print('##', minimizer_record_file1, discard_mini_mer[:100], 'tt_3mer_num', int(tt_3mer_num / (10 ** 3)), max(min(200, int(tt_3mer_num / (10 ** 3))), 20), 'len(minimizer_list)', len(minimizer_list))
	discard_mini_mer_set = set([ref_mini_3mer for ref_mini_3mer,mer_num in discard_mini_mer if mer_num > len(minimizer_list) / 100])
	#print(minimizer_record_file1, discard_mini_mer_set, len(list(discard_mini_mer_set)))
	read_chunk = []
	chunk_number = 0
	output_clust_data = []
	keep_reads = []
	total_seq_num = len(ref_name_list)
	line_num = 0
	#print('discard 3mer finished')
	#Math read with mini-3mer index
	for data_line in minimizer_list:
		start_time = time()
		line_num += 1
		#print('line_num', line_num, line_num in clustered_reads)
		output_data_line = []
		if line_num in clustered_reads:
			continue
		reads_in_cluster = []
		#reads_in_cluster.append((line_num, minimizer_list[line_num - 1], seq_list[line_num - 1], 10 ** 4))
		reads_in_cluster.append((line_num, minimizer_list[line_num - 1],  10 ** 4))
		query_name, query_namev, query_minimizers_lst, query_minimizers_lstv = decord_minimizers(data_line)
		match_count = Counter()
		query_minimizers_lst_len = len(query_minimizers_lst)

		#Match the foward sequence
		counter_time = time()
		mini_3mer_lst = []
		discard_mer_num = 0
		
		
		for m in range(query_minimizers_lst_len - 2):
			minimizer1, minimizer2, minimizer3 = query_minimizers_lst[m], query_minimizers_lst[m + 1], query_minimizers_lst[m + 2]
			query_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
			mini_3mer_lst.append(query_mini_3mer)
		for query_mini_3mer in list(set(mini_3mer_lst)):
			if query_mini_3mer in discard_mini_mer_set:
				discard_mer_num += 1
				continue
			if query_mini_3mer in mini_mer_index:
				match_count.update(list(set(mini_mer_index[query_mini_3mer])))
				#print(len(list(set(mini_mer_index[query_mini_3mer]))))

		query_seq = seq_dic[line_num]
		query_7mer_set = set([query_seq[i : i + 7] for i in range(len(query_seq) - 7 + 1)])

		match_count_lst = list(match_count.most_common())
		#print('match_info', len(list(set(mini_3mer_lst))), discard_mer_num)
		match_num_cutoff = 1 + int(0.1 * (len(list(set(mini_3mer_lst))) - discard_mer_num))
		for match_result in match_count_lst:
			ref_record_num, ref_match_num = match_result
			if ref_record_num == line_num:
				continue
			if ref_match_num < match_num_cutoff:
				break
			ref_seq = seq_dic[ref_record_num]
			dif_transcript_ratio = confirm_not_different_transcript(ref_seq, query_7mer_set)
			if dif_transcript_ratio > 0.75:
				#print('dif_transcript_ratio, ref_match_num', dif_transcript_ratio, ref_match_num, match_num_cutoff, ref_name_list[line_num], ref_name_list[ref_record_num])
				#if int(ref_name_list[line_num] / 100000) != int(ref_name_list[ref_record_num] / 100000):
				#	print(seq_dic[line_num], seq_dic[ref_record_num])
				reads_in_cluster.append((ref_record_num, minimizer_list[ref_record_num - 1], ref_match_num))
		#Match the reverse sequence
		match_count = Counter()
		query_minimizers_lstv_len = len(query_minimizers_lstv)
		counter_time = time()
		mini_3mer_lst = []
		discard_mer_num = 0
		for m in range(query_minimizers_lstv_len - 2):
			minimizer1, minimizer2, minimizer3 = query_minimizers_lstv[m], query_minimizers_lstv[m + 1], query_minimizers_lstv[m + 2]
			query_mini_3mer = generate_3mer(minimizer1, minimizer2, minimizer3)
			mini_3mer_lst.append(query_mini_3mer)
		for query_mini_3mer in list(set(mini_3mer_lst)):
			if query_mini_3mer in discard_mini_mer_set:
				discard_mer_num += 1
				continue
			if query_mini_3mer in mini_mer_index:
				match_count.update(list(set(mini_mer_index[query_mini_3mer])))
		match_count_lstv = list(match_count.most_common())
		match_num_cutoffv = 1 + int(0.1 * (len(list(set(mini_3mer_lst))) - discard_mer_num))
		#print('match_info', len(list(set(mini_3mer_lst))), discard_mer_num)

		query_seqv = seq_dic[-line_num]
		query_7mer_set = set([query_seq[i : i + 7] for i in range(len(query_seq) - 7 + 1)])

		for match_resultv in match_count_lstv:
			ref_record_numv, ref_match_numv = match_resultv
			if ref_match_numv < match_num_cutoffv:
				break
			ref_seqv = seq_dic[ref_record_numv]
			dif_transcript_ratiov = confirm_not_different_transcript(ref_seqv, query_7mer_set)
			if dif_transcript_ratiov > 0.75:
				#print('dif_transcript_ratiov, ref_match_num', dif_transcript_ratiov, ref_match_numv, match_num_cutoffv, len(list(set(mini_3mer_lst))), ref_name_list[line_num], ref_name_list[ref_record_numv])
				reads_in_cluster.append((-ref_record_numv, minimizer_list[ref_record_numv - 1], ref_match_numv))

		#start_time = time()
		#Find the best read in this cluster and record
		best_read_in_cluster, reads_in_cluster_refine = find_best_reads(reads_in_cluster)
		#print('refine time', time() - start_time, len(reads_in_cluster_refine))
		
		output_data_line.append('{}\t'.format(int(best_read_in_cluster/abs(best_read_in_cluster) * ref_name_list[abs(best_read_in_cluster)])))
		for read_name_in_cluster, mer_num in reads_in_cluster_refine:
			output_data_line.append('{}|{}\t'.format(int(read_name_in_cluster/abs(read_name_in_cluster) * ref_name_list[abs(read_name_in_cluster)]), mer_num))
			clustered_reads.add(abs(read_name_in_cluster))
		clustered_reads.add(abs(best_read_in_cluster))
		keep_reads.append(abs(best_read_in_cluster))
		if output_data_line != []:
			output_clust_data.append(''.join((''.join(output_data_line), '\n')))
		#print(''.join((''.join(output_data_line), '\n')), line_num)
	keep_reads = list(set(keep_reads))
	keep_reads.sort()
	keep_data = [minimizer_list[i - 1] for i in keep_reads]
	seq_keep_data = [seq_lst[i - 1] for i in keep_reads]
	keep_data_file = open('{}/minimizer_record_sort/{}'.format(outdir, minimizer_record_file1), 'w')
	keep_data_file.writelines(keep_data)
	keep_data_file.close()
	seq_keep_data_file = open('{}/seq_record_sort/{}'.format(outdir, minimizer_record_file1), 'w')
	seq_keep_data_file.writelines(seq_keep_data)
	seq_keep_data_file.close()
	clust_record_file = open('{}/match_record/{}-{}_{}'.format(outdir, minimizer_record_file1, minimizer_record_file1, r_n), 'w')
	clust_record_file.writelines(output_clust_data)
	clust_record_file.close()
	print("\033[31m[Cluster sequence]\033[0m", 'finished self compare {}, round {}'.format(minimizer_record_file1, r_n), 'time_used: %.3f s'%( time() - compare_start_time))



def cluster_reads_1st(args, outdir):
	threads = args.t
	os.system("rm -rf {}/match_record".format(outdir))
	os.mkdir('{}/match_record'.format(outdir))

	print("\033[31m**************\033[0m", 'Clustering using minimizer-3mer', "\033[31m**************\033[0m")

	start_time = time()
	minimizer_record_file_lst = [int(w) for w in os.listdir('{}/minimizer_record_sort'.format(outdir))]

	minimizer_record_file_lst.sort()

	clust_p = mp.Pool(processes = threads)
	for minimizer_record_file1 in minimizer_record_file_lst:
		clust_p.apply_async(compare_file_self, (minimizer_record_file1, 1, outdir, ))
	clust_p.close()
	clust_p.join()

	clust_p = mp.Pool(processes = threads)
	for minimizer_record_file1 in minimizer_record_file_lst:
		clust_p.apply_async(compare_file_self, (minimizer_record_file1, 2, outdir, ))
	clust_p.close()
	clust_p.join()

	clust_p = mp.Pool(processes = threads)
	for minimizer_record_file1 in minimizer_record_file_lst:
		num_len_line1 = open('{}/minimizer_record_len/{}'.format(outdir, minimizer_record_file1)).readlines()[0].split()
		file1_minimizer_num_h, file1_minimizer_num_l = int(num_len_line1[0]), int(num_len_line1[1])
		for minimizer_record_file2 in minimizer_record_file_lst:
			num_len_line2 = open('{}/minimizer_record_len/{}'.format(outdir, minimizer_record_file2)).readlines()[0].split()
			file2_minimizer_num_h, file2_minimizer_num_l = int(num_len_line2[0]), int(num_len_line2[1])
			if int(minimizer_record_file2) <= int(minimizer_record_file1) or int(minimizer_record_file2) >= int(minimizer_record_file1) + max(10, len(minimizer_record_file_lst) / 50):
				continue
			if int(minimizer_record_file2) == int(minimizer_record_file1) + 1 or file2_minimizer_num_l / file1_minimizer_num_h > 0.7:
				#print(minimizer_record_file1 , minimizer_record_file2)
				clust_p.apply_async(compare_file2file, (minimizer_record_file1, minimizer_record_file2, outdir, ))
	clust_p.close()
	clust_p.join()
	print("\033[31m[Cluster sequence finished]\033[0m", 'finished read cluster, time used: {}  s'.format(time() - start_time))

	

