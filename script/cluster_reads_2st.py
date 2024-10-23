#!/data/software/conda/bin/python
import multiprocessing as mp
import os,sys
from time import time
import math
from collections import deque, Counter
import gzip
import parasail
import re



def rev_string(string):
	return ''.join([string[len(string)-n-1] for n in range(0,len(string))])

def rev_com_seq(sequence):
	nuc_rev={'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a','N':'N','X':'X','n':'n','Y':'R','R':'Y','K':'M','M':'K','S':'S','W':'W','B':'V','V':'B','H':'D','D':'H','y':'r','r':'y','k':'m','m':'k','s':'s','w':'w','b':'v','v':'b','h':'d','d':'h', '0': '3', '3': '0', '2': '1', '1': '2', '4': '4'} 
	return ''.join([nuc_rev[sequence[len(sequence)-n-1]] for n in range(0,len(sequence))])
 
def minimizer_poisson_mean(qual):
	return (sum((1 - D[char_] for char_ in qual))/len(qual))

def get_minimizers_seq(sequence, k, w):
	minimizers = deque([(sequence[i: i + k], i) for i in range(w - k + 1) if i + k <= len(sequence)])
	if len(minimizers) == 0:
		return []
	minimizer, minimizers_hash = min(minimizers), [min(minimizers)]
	for i in range(w - k + 1,len(sequence) - k):
		discard_kmer, new_kmer = minimizers.popleft(), (sequence[i: i + k], i)
		minimizers.append(new_kmer)
		if discard_kmer == minimizer:
			minimizer = min(minimizers)
			minimizers_hash.append(minimizer)
		elif new_kmer < minimizer:
			minimizer = new_kmer
			minimizers_hash.append(minimizer)      
	#print(minimizers_hash) 
	return [int(minimizer) for minimizer, lo in minimizers_hash]

def cigar_seqs(sequence1, sequence2):
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
	return cigar_result


def seq_similarity_coverage(sequence1, sequence2):
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
	for cig in cigar_result:
		if cig[1] == '=':
			nt_match += cig[0]
		else:
			nt_mismatch += cig[0]
		if cig[1] == '=' or cig[1] == 'I':
			sequence2_in_align += cig[0]
	if nt_match+nt_mismatch == 0:
		return 0, 0
	#if nt_match < len(sequence1)/3 or nt_match < len(sequence2)/3:
		#return 0
	if nt_match < 20:
		return 0, 0
	#print(nt_match/(nt_match+nt_mismatch))
	return nt_match/(nt_match+nt_mismatch), sequence2_in_align / len(sequence1)


def SW_match(read_minimizer_info):
	line_num, ref_seq, ref_seq_rev, ref_minimizers, ref_minimizers_rev = read_minimizer_info
	read_counter = Counter([])
	match_num_cutoff = 1 + int(0.1 * (len(list(set(ref_minimizers)))))
	related_reads = []
	for minimizer in ref_minimizers:
		if minimizer in MINIMIZER_INDEX:
			read_counter.update(MINIMIZER_INDEX[minimizer])
	read_best_support = list(read_counter.most_common(5))
	#print(line_num, read_best_support[:10])
	for query_line_num, match_num in read_best_support:
		if match_num < match_num_cutoff:
			break
		if line_num == query_line_num:
			continue
		query_seq = MINIMIZER_DATA_CHUNK[query_line_num - 1 - MAP_CHUNK_ORDER][1]
		if len(query_seq) < len(ref_seq):
			continue
		match_identity, match_coverage = seq_similarity_coverage(ref_seq, query_seq)
		#if match_identity > 0.75 and match_coverage > 0.5:
		#print(match_identity, match_coverage)
		#print(cigar_seqs(ref_seq, query_seq))
		if match_identity > 0.85 and match_coverage > 0.85:
			related_reads.append((line_num, query_line_num))
	
	read_counter_rev = Counter([])
	match_num_cutoff_rev = 1 + int(0.1 * (len(list(set(ref_minimizers_rev)))))
	for minimizer in ref_minimizers_rev:
		if minimizer in MINIMIZER_INDEX:
			read_counter_rev.update(MINIMIZER_INDEX[minimizer])
	read_best_support_rev = list(read_counter_rev.most_common(5))
	#print(-line_num, read_best_support_rev[:10])
	for query_line_num, match_num in read_best_support_rev:
		if match_num < match_num_cutoff_rev:
			break
		if line_num == query_line_num:
			continue
		query_seq = MINIMIZER_DATA_CHUNK[query_line_num - 1 - MAP_CHUNK_ORDER][1]
		if len(query_seq) < len(ref_seq_rev):
			continue
		match_identity, match_coverage = seq_similarity_coverage(ref_seq_rev, query_seq)
		#print(match_identity, match_coverage)
		#print(cigar_seqs(ref_seq_rev, query_seq))
		if match_identity > 0.85 and match_coverage > 0.85:
			#print(query_seq, ref_seq_rev)
			#print(cigar_seqs(ref_seq_rev, query_seq))
			related_reads.append((-line_num, query_line_num))
	return related_reads

def confirm_not_different_transcript(query_seq, ref_7mer_lst):
	query_7mer_set = set([query_seq[i : i + 7] for i in range(len(query_seq) - 7 + 1)])
	ref_7mer_in_num = 0
	for kmer in ref_7mer_lst:
		if kmer in query_7mer_set:
			ref_7mer_in_num += 1
	return ref_7mer_in_num / (len(ref_7mer_lst) + 1)

def K_match(read_minimizer_info):
	line_num, ref_seq, ref_seq_rev, ref_minimizers, ref_minimizers_rev = read_minimizer_info
	read_counter = Counter([])
	match_num_cutoff = 1 + int(0.1 * (len(list(set(ref_minimizers)))))
	related_reads = []
	for minimizer in ref_minimizers:
		if minimizer in MINIMIZER_INDEX:
			read_counter.update(MINIMIZER_INDEX[minimizer])
	read_best_support = list(read_counter.most_common())
	ref_7mer_lst = [ref_seq[i : i + 7] for i in range(len(ref_seq) - 7 + 1)]
	#print(line_num, read_best_support[:10])
	for query_line_num, match_num in read_best_support:
		if match_num < match_num_cutoff:
			break
		if line_num == query_line_num:
			continue
		query_seq = MINIMIZER_DATA_CHUNK[query_line_num - 1 - MAP_CHUNK_ORDER][1]
		if len(query_seq) < len(ref_seq):
			continue
		#match_identity, match_coverage = seq_similarity_coverage(ref_seq, query_seq)
		#if match_identity > 0.75 and match_coverage > 0.5:
		#print(match_identity, match_coverage)
		#print(cigar_seqs(ref_seq, query_seq))
		#if match_identity > 0.85 and match_coverage > 0.85:
		#	related_reads.append((line_num, query_line_num))
		dif_transcript_ratio = confirm_not_different_transcript(query_seq, ref_7mer_lst)
		#print(match_identity, match_coverage, dif_transcript_ratio)
		if dif_transcript_ratio > 0.8:
			related_reads.append((line_num, query_line_num))

	read_counter_rev = Counter([])
	match_num_cutoff_rev = 1 + int(0.1 * (len(list(set(ref_minimizers_rev)))))
	for minimizer in ref_minimizers_rev:
		if minimizer in MINIMIZER_INDEX:
			read_counter_rev.update(MINIMIZER_INDEX[minimizer])
	read_best_support_rev = list(read_counter_rev.most_common())
	ref_7mer_lst_rev = [ref_seq_rev[i : i + 7] for i in range(len(ref_seq_rev) - 7 + 1)]
	#print(-line_num, read_best_support_rev[:10])
	for query_line_num, match_num in read_best_support_rev:
		if match_num < match_num_cutoff_rev:
			break
		if line_num == query_line_num:
			continue
		query_seq = MINIMIZER_DATA_CHUNK[query_line_num - 1 - MAP_CHUNK_ORDER][1]
		if len(query_seq) < len(ref_seq_rev):
			continue

		#match_identity, match_coverage = seq_similarity_coverage(ref_seq_rev, query_seq)
		#print(match_identity, match_coverage)
		#print(cigar_seqs(ref_seq_rev, query_seq))
		#if match_identity > 0.85 and match_coverage > 0.85:
			#print(query_seq, ref_seq_rev)
			#print(cigar_seqs(ref_seq_rev, query_seq))
			#related_reads.append((-line_num, query_line_num))
		dif_transcript_ratio = confirm_not_different_transcript(query_seq, ref_7mer_lst_rev)
		#print(match_identity, match_coverage, dif_transcript_ratio)
		if dif_transcript_ratio > 0.8:
			related_reads.append((-line_num, query_line_num))
	return related_reads

def initializer_SW_match(minimizer_data_chunk, minimizer_index, map_chunk_order):
	global MINIMIZER_DATA_CHUNK, MINIMIZER_INDEX, MAP_CHUNK_ORDER
	MINIMIZER_DATA_CHUNK, MINIMIZER_INDEX, MAP_CHUNK_ORDER = minimizer_data_chunk, minimizer_index, map_chunk_order

def initializer_K_match(minimizer_data_chunk, minimizer_index, map_chunk_order):
	global MINIMIZER_DATA_CHUNK, MINIMIZER_INDEX, MAP_CHUNK_ORDER
	MINIMIZER_DATA_CHUNK, MINIMIZER_INDEX, MAP_CHUNK_ORDER = minimizer_data_chunk, minimizer_index, map_chunk_order

def generate_minimizer_data(read_info):
	line_num, ref_seq = read_info
	ref_seq_rev = rev_com_seq(ref_seq)
	ref_minimizers = get_minimizers_seq(ref_seq, k, w)
	ref_minimizers_rev = get_minimizers_seq(ref_seq_rev, k, w)
	return line_num, ref_seq, ref_seq_rev, ref_minimizers, ref_minimizers_rev
	

def get_new_cluster_id(cluster_num_set):
	for i in range(10**100):
		if i not in cluster_num_set:
			return i
			break

def refine_clusters(refined_relations_uniform):

	reads_clustered = {}
	cluster_lst_dic = {}
	cluster_num_set = set([])
	for query_read, ref_read in refined_relations_uniform:
		bst_read, match_read = int(ref_read), int(query_read)
		bst_read_abs, match_read_abs = abs(int(ref_read)), abs(int(query_read))

		if bst_read_abs not in reads_clustered and match_read_abs not in reads_clustered:
			#print('Both not in')
			new_cluster_id = get_new_cluster_id(cluster_num_set)
			cluster_num_set.add(new_cluster_id)
			reads_clustered[bst_read_abs] = (new_cluster_id, int(bst_read_abs / bst_read))
			reads_clustered[match_read_abs] = (new_cluster_id, int(match_read_abs / match_read))
			cluster_lst_dic[new_cluster_id] = [bst_read, match_read]
			continue
		elif bst_read_abs not in reads_clustered and match_read_abs in reads_clustered:
			#print('Bst not in')
			match_cluster_id, match_cluster_fr = reads_clustered[match_read_abs]
			fow_bac_mark = int((bst_read_abs/ bst_read) * (match_read_abs/ match_read) * match_cluster_fr)
			#print('fow_bac_mark', fow_bac_mark)
			reads_clustered[bst_read_abs] = (match_cluster_id, fow_bac_mark)
			cluster_lst_dic[reads_clustered[match_read_abs][0]].append(bst_read_abs * fow_bac_mark)
		elif bst_read_abs in reads_clustered and match_read_abs not in reads_clustered:
			#print('Match not in')
			bst_cluster_id, bst_cluster_fr = reads_clustered[bst_read_abs]
			fow_bac_mark = int((bst_read_abs/ bst_read) * (match_read_abs/ match_read) * bst_cluster_fr)
			#print('fow_bac_mark', fow_bac_mark)
			reads_clustered[match_read_abs] = (bst_cluster_id, fow_bac_mark)
			cluster_lst_dic[reads_clustered[bst_read_abs][0]].append(match_read_abs * fow_bac_mark)
		elif  bst_read_abs in reads_clustered and match_read_abs in reads_clustered:
			bst_cluster_id, bst_cluster_fr = reads_clustered[bst_read_abs]
			match_cluster_id, match_cluster_fr = reads_clustered[match_read_abs]
			if bst_cluster_id == match_cluster_id:
				#print('Both in, same cluster')
				continue
			if bst_cluster_id != match_cluster_id:
				#print('Both in, dif cluster')
				fow_bac_mark = int((bst_read_abs/ bst_read) * (match_read_abs/ match_read) * bst_cluster_fr * match_cluster_fr)
				#print('fow_bac_mark', fow_bac_mark)
				cluster_lst_dic[bst_cluster_id] += [int(abs(w) * fow_bac_mark) for w in cluster_lst_dic[match_cluster_id]]
				for read_name in cluster_lst_dic[match_cluster_id]:
					#reads_clustered[abs(read_name)] = (bst_cluster_id, reads_clustered[abs(match_read)][1] * fow_bac_mark)
					reads_clustered[abs(read_name)] = (bst_cluster_id, fow_bac_mark)
				del cluster_lst_dic[match_cluster_id]
				cluster_num_set.remove(match_cluster_id)
			continue
		#print(cluster_lst_dic[reads_clustered[bst_read_abs][0]])
	return cluster_lst_dic


def cluster_reads_2st(args, outdir):
	start_time = time()
	start_time_t = time()
	print("\033[31m**************\033[0m", 'Refine clusters', "\033[31m**************\033[0m")
	global k, w
	k = args.k + 2
	w = args.w
	threads = args.t
	clusters_1st = {}
	ref_seq_dic = {}
	ref_minimizer_lst = []
	line_num = 0
	ref_seq_list = []
	with open ('{}/ref_seqs_1st.txt'.format(outdir)) as cluster_1st_data:
		for line in cluster_1st_data:
			line_num += 1
			line_splt = line.split()
			ref_seq = line_splt[1]
			clusters_1st[line_num] = line_splt[2:]
			ref_seq_list.append((line_num, ref_seq))
	print("\033[31m[Read cluster results]\033[0m", 'finished get representative read, time used: {}  s'.format(time() - start_time))

	start_time = time()
	p = mp.Pool(processes = threads)
	minimizer_data = p.map(generate_minimizer_data, ref_seq_list)
	p.close()
	p.join()
	minimizer_data.sort(key = lambda x:x[0])
	print("\033[31m[Generate Kmers]\033[0m", 'finished get Kmers, time used: {}  s'.format(time() - start_time))

	refined_relations_uniform = []
	map_chunk_size = 10 ** 5
	for chunk_n in range(int((len(minimizer_data) - 1) /map_chunk_size  + 1)):
		start_time = time()
		minimizer_index = {}
		#print(int((len(minimizer_data) - 1) / 10 ** 5 + 1), len(minimizer_data))
		minimizer_data_chunk = minimizer_data[chunk_n * (map_chunk_size) : min((chunk_n + 1) * (map_chunk_size), len(minimizer_data))]
		for line_num, ref_seq, ref_seq_rev, ref_minimizers, ref_minimizers_rev in minimizer_data_chunk:
			for minimizer in ref_minimizers:
				if minimizer not in minimizer_index:
					minimizer_index[minimizer] = []
				minimizer_index[minimizer].append(line_num)
		discard_minimizer = [(0, 0) for i in range(500)]
		dele_ref_minimizers = []
		tt_minimizer_num = 0
		for i, minimizer in enumerate(minimizer_index):
			minimizer_num = len(minimizer_index[minimizer])
			#if minimizer_num == 1:
			#	dele_ref_minimizers.append(minimizer)
			#	continue
			tt_minimizer_num += 1
			if minimizer_num > discard_minimizer[-1][1]:
				discard_minimizer[-1] = (minimizer, minimizer_num)
				discard_minimizer.sort(key = lambda x:x[1], reverse = True)
		#for del_minimizer in dele_ref_minimizers:
		#	del minimizer_index[del_minimizer]
		discard_minimizer = discard_minimizer[: min([max(min(200, int(tt_minimizer_num / (10 ** 3))), 20), int(tt_minimizer_num / (10))])]
		#print(discard_minimizer)
		discard_minimizer_lst = list(set([ref_minimizer for ref_minimizer,minimizer_read_num in discard_minimizer if minimizer_read_num > len(minimizer_data_chunk) / 100]))
		for del_minimizer in discard_minimizer_lst:
			del minimizer_index[del_minimizer]

		print("\033[31m[Generate K-mers]\033[0m", 'finished  construct {}-mer index in chunk {} time used: {}  s'.format(k, chunk_n + 1, time() - start_time))

		######################single thread for debug##############
		#initializer_SW_match(minimizer_data, minimizer_index)
		#for read_info in minimizer_data:
		#	SW_match(read_info)
		###########################################################
		start_time = time()
		map_chunk_order = chunk_n * map_chunk_size
		p = mp.Pool(processes = threads, initializer = initializer_K_match(minimizer_data_chunk, minimizer_index, map_chunk_order))
		refined_relations = p.map(SW_match, minimizer_data)
		p.close()
		p.join()
	
		for rela_data in refined_relations:
			#if rela_data != []:
			#print(rela_data)
			for rela in rela_data:
				#if rela[0] < 0 or rela[1] < 0:
					#print(rela)
				refined_relations_uniform.append(rela)
		print("\033[31m[K-mer seed map]\033[0m", 'finished {}-mer mapping in chunk {}, got {} relations, time used: {}  s'.format(k, chunk_n + 1, len(refined_relations_uniform),  time() - start_time))


	start_time = time()
	#print(len(refined_relations_uniform))
	refined_cluster = []
	refined_cluster_set = set([])
	cluster_lst_dic = refine_clusters(refined_relations_uniform)
	#print('cluster_lst_dic', cluster_lst_dic)	
	for _, refine_i in enumerate(cluster_lst_dic):
		refine_line_lst = cluster_lst_dic[refine_i]
		#print([minimizer_data[abs(t) - 1][1] for t in cluster_lst_dic[i]])
		refine_line_lst.sort(key = lambda x:abs(x))
		#print(refine_line_lst)
		best_read = ''
		best_read_seq = ''
		reads_in_refined_cluster = []
		for line_num in refine_line_lst:
			line_fr = int(abs(line_num) / line_num)
			line_num_abs = abs(line_num)
			if best_read == '':
				best_read = int(clusters_1st[line_num_abs][0].split('|')[0]) * line_fr
				line_num, ref_seq, ref_seq_rev, ref_minimizers, ref_minimizers_rev = minimizer_data[line_num_abs - 1]
				if line_fr == 1:
					best_read_seq = ref_seq
				if line_fr == -1:
					best_read_seq = ref_seq_rev
			#print(clusters_1st[line_num_abs])
			reads_in_refined_cluster += [int(read_id.split('|')[0]) * line_fr for read_id in clusters_1st[line_num_abs]]
			refined_cluster_set.add(line_num_abs)
		#print(best_read, best_read_seq, reads_in_refined_cluster)
		refined_cluster.append((best_read, best_read_seq, reads_in_refined_cluster))
	
	for line_num, ref_seq, ref_seq_rev, ref_minimizers, ref_minimizers_rev in minimizer_data:
		if line_num in refined_cluster_set:
			continue
		line_num_abs = abs(line_num)
		best_read = int(clusters_1st[line_num_abs][0].split('|')[0])
		best_read_seq = ref_seq
		reads_in_refined_cluster = [int(read_id.split('|')[0]) for read_id in clusters_1st[line_num_abs]]
		refined_cluster.append((best_read, best_read_seq, reads_in_refined_cluster))
	refined_cluster.sort(key = lambda x:len(x[2]), reverse = True)
	output_file = open('{}/ref_seqs_2st.txt'.format(outdir), 'w')
	for best_read, read_seq, read_lst in refined_cluster:
		output_line = ''.join(('{}\t{}\t'.format(best_read, read_seq), ''.join([str(read_id) + '\t' for read_id in read_lst]), '\n'))
		output_file.writelines(output_line)
	output_file.close()
	
	
	#print('refined cluster number', len(refined_cluster))


	print("\033[31m[Refine clusters]\033[0m", 'finished refine cluster, time used: {}  s'.format(time() - start_time_t))

	

