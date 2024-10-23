#!/data/software/conda/bin/python


import networkx as nx
import multiprocessing as mp
import os, sys
import random
from collections import deque, Counter
import parasail
import re
from time import time




def find_best_seq(seq_lst):
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
	seq_lst_sorted = [(seq_name, seq_score) for seq_name, seq, seq_score in seq_lst_scored]
	return best_seq_name, best_seq, seq_lst_sorted

def refine_read_lst(read_lst_lst):
	tt_read_lst = []
	for read_lst in read_lst_lst:
		read_abs_lst = [abs(read_name) for read_name in read_lst]
		tt_read_lst += read_abs_lst
		#if len(read_abs_lst) !=  len(list(set(read_abs_lst))):
		#	print('!!!!!!!!!!!!!!!!!!!!!!repeat read name in cluster')
	#print(len(tt_read_lst), len(list(set(tt_read_lst))))
	return read_lst_lst

def get_gene_list_from_graph(graph_lst):
	gene_dic = {}
	while graph_lst != []:
		remove_edges = []
		for edge_ in graph_lst:
			read1, read2, connect_weight = edge_
			if gene_dic == {}:
				gene_dic[read1] = 1
				gene_dic[read2] = gene_dic[read1] * connect_weight
				remove_edges.append(edge_)
				continue
			else:
				if read1 not in gene_dic and read2 not in gene_dic:
					continue
				if read1 in gene_dic and read2 in gene_dic:
					#print('?' * 10)
					remove_edges.append(edge_)
				if read1 in gene_dic and read2 not in gene_dic:
					gene_dic[read2] = gene_dic[read1] * connect_weight
					remove_edges.append(edge_)
				if read1 not in gene_dic and read2 in gene_dic:
					gene_dic[read1] = gene_dic[read2] * connect_weight
					remove_edges.append(edge_)
		if len(remove_edges) == 0:
			break
		remove_edges_set = set(remove_edges)
		graph_lst = [g for g in graph_lst if g not in remove_edges_set]
	gene_lst = [gene * gene_dic[gene] for i, gene in enumerate(gene_dic)]
	return gene_lst	

def extract_reads_from_files(read_list, chunk_read_num):
	reads_file_dic = {}
	reads_in_chunk = [(read_list[i], i) for i in range(len(read_list))]
	for seq_num_raw, order_t in reads_in_chunk:
		seq_num = abs(seq_num_raw)
		file_name, record_num = int((seq_num -1) / chunk_read_num), (seq_num - 1) % chunk_read_num + 1
		if file_name not in reads_file_dic:
			reads_file_dic[file_name] = []
		reads_file_dic[file_name].append((record_num, order_t, seq_num_raw))
	tt = 0
	reads_file_lst = []
	for i,file_name in enumerate(reads_file_dic):
		reads_file_lst.append((file_name, reads_file_dic[file_name]))
	random.shuffle(reads_file_lst)

	max_minimizer_num = 0
	min_minimizer_num = 10 ** 10
	seq_chunk_data = []
	for file_name,read_record_lst in reads_file_lst:
		#print('scanned read file {}'.format(file_name))
		read_record_lst.sort(key = lambda x:x[0])
		tt = 0
		read_record_deque = deque(read_record_lst)
		if len(read_record_deque) != 0:
			fst_read_num, fst_read_order, fst_seq_num = read_record_deque.popleft()
		else:
			break
		line_num = 0
		with open('{}/seq_record/{}'.format(OUTDIR, file_name)) as seq_record_file:
			for seq_record_line in seq_record_file:
				line_num +=1
				if line_num == fst_read_num:
					compressed_seq, compressed_seqv = seq_record_line.split('\t')[1], seq_record_line.split('\t')[3]
					if fst_seq_num > 0:
						seq_chunk_data.append((compressed_seq, fst_read_order, fst_seq_num))
					if fst_seq_num < 0:
						seq_chunk_data.append((compressed_seqv, fst_read_order, fst_seq_num))
					if len(read_record_deque) == 0:
						break
					else:
						fst_read_num, fst_read_order, fst_seq_num = read_record_deque.popleft()
	
	#seq_chunk_data.sort(key = lambda x:x[1])
	seq_chunk_output_data = {seq_num: seq for seq, seq_order, seq_num in seq_chunk_data}
	return seq_chunk_output_data

def get_new_cluster_id(cluster_num_set):
	for i in range(10**100):
		if i not in cluster_num_set:
			return i
			break

def combine_file2file_results(file_name):
	#print('###', file_name)
	reads_clustered = {}
	cluster_lst_dic = {}
	cluster_num_set = set([])
	match_in_chunk = {}
	file_mark = ''.join(('-', str(file_name)))
	fl = len(file_mark)
	file_need = []
	for match_file in match_file_lst:
		if '_' in match_file:
			continue
		if match_file[-fl:] == file_mark:
			file_need.append((match_file, int(match_file.split('-')[0])))
		file_need.sort(key = lambda x:x[1], reverse = True)
	for match_file,_ in file_need:
		#print(file_name, match_file)
		with open('{}/match_record/{}'.format(OUTDIR, match_file)) as match_file_data:
			for match_line in match_file_data:
				#print('##########match line', match_line.split())
				for match_term in match_line.split():
					#print('mtch_term', match_term)
					ref_read, query_read, match_num = match_term.split('|')
					if int(match_num) <= 1:
						continue
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
					
	#cluster_lst = [cluster_lst_dic[cluster_id] for i,cluster_id in enumerate(cluster_lst_dic)]
	#cluster_lst.sort(key = lambda x:len(x), reverse = True)
	#print(len(cluster_lst))
	return file_name, reads_clustered, cluster_num_set, cluster_lst_dic

def add_self_cluster(chunk_cluster):
	file_name, reads_clustered, cluster_num_set, cluster_lst_dic = chunk_cluster
	for rnd in [1, 2]:
		fst_results = open('{}/match_record/{}-{}_{}'.format(OUTDIR, file_name, file_name, rnd)).readlines()

		for fst_line in fst_results:
			fst_line_splt = fst_line.split()
			bst_read = int(fst_line_splt[0])
			bst_read_abs = abs(bst_read)
			if len(fst_line_splt) == 1:
				continue
			for match_info in fst_line_splt[1:]:
				#print(match_info)
				match_read, match_num = int(match_info.split('|')[0]), float(match_info.split('|')[1])
				match_read_abs = abs(match_read)
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
	#cluster_lst = [cluster_lst_dic[cluster_id] for i,cluster_id in enumerate(cluster_lst_dic)]
	#cluster_lst.sort(key = lambda x:len(x), reverse = True)
	#print(file_name, len(cluster_lst), cluster_lst[-10 : ])
	return file_name, reads_clustered, cluster_num_set, cluster_lst_dic

def combine_2_chunk_cluster_lst(clust_chunk_tpl):
	if len(clust_chunk_tpl) == 1:
		return clust_chunk_tpl[0]
	cluster_data1, cluster_data2 = clust_chunk_tpl
	file_name, reads_clustered, cluster_num_set, cluster_lst_dic = cluster_data1
	_, _, _, cluster_to_combine = cluster_data2
	for i, cluster_id in enumerate(cluster_to_combine):
		read_lst_to_combine = cluster_to_combine[cluster_id]
		bst_read = int(read_lst_to_combine[0])
		bst_read_abs = abs(bst_read)
		for match_read in read_lst_to_combine[1:]:
			match_read_abs = abs(match_read)
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
	return file_name, reads_clustered, cluster_num_set, cluster_lst_dic

def cluster_reads_1st_combine(args, outdir):
	start_time = time()
	print("\033[31m**************\033[0m", 'Connecting reads in graph', "\033[31m**************\033[0m")
	global OUTDIR
	OUTDIR = outdir
	threads = args.t
	read_chunk_size = args.c
	clust_99_in = {}
	clust_99_out = {}

	file_lst = [int(f) for f in os.listdir('{}/minimizer_record_len'.format(outdir))]
	file_lst.sort(reverse = True)
	
	file_lst2 = [int(f) for f in os.listdir('{}/minimizer_record_len'.format(outdir))]
	file_lst2.sort(reverse = True)

	#print('file_lst', file_lst)
	global match_file_lst
	match_file_lst = os.listdir('{}/match_record'.format(outdir))
	#for file_name in file_lst:	
	p = mp.Pool(processes = threads)
	chunk_cluster_lst = p.map(combine_file2file_results, file_lst)
	p.close()
	p.join

	p = mp.Pool(processes = threads)
	chunk_cluster_lst = p.map(add_self_cluster, chunk_cluster_lst)
	p.close()
	p.join

	chunk_cluster_lst.sort(key = lambda x:x[0])

	NN_G = nx.Graph()

	path_set = set([])
	for file_name, reads_clustered, cluster_num_set, cluster_lst_dic in chunk_cluster_lst:
		for i, cluster_id in enumerate(cluster_lst_dic):
			read_lst_to_combine = cluster_lst_dic[cluster_id]
			bst_read = int(read_lst_to_combine[0])
			bst_read_abs = abs(bst_read)
			for match_read in read_lst_to_combine[1:]:
				match_read_abs = abs(match_read)
				path_weight = int((bst_read_abs * match_read_abs) / (bst_read * match_read))
				path = (bst_read_abs, match_read_abs, path_weight)
				path2 = (match_read_abs,bst_read_abs, path_weight)
				if path in path_set or path2 in path_set:
					continue
				
				path_set.add(path)
				NN_G.add_edge( bst_read_abs, match_read_abs, weight = path_weight)

	sub_nets = nx.connected_components(NN_G)
	sub_graphs = list([NN_G.subgraph(c) for c in nx.connected_components(NN_G)])
	print("\033[31m[Graph finished]\033[0m", 'consctruced {} connected subgraphs, time used: {}  s'.format(len(sub_graphs), time() - start_time))
	start_time = time()
	print("\033[31m**************\033[0m", 'Get best read in each graph', "\033[31m**************\033[0m")
	z = 1
	ref_seqs = []
	ref_seqs_cluster = []
	threads = 10
	read_lst_lst = []
	tt_read_lst = []

	for w in sub_graphs:
		weight_lst = list(w.edges(data = 'weight'))
		#if len(weight_lst) >= 3:
		read_lst = get_gene_list_from_graph(weight_lst)
		read_lst_lst.append(read_lst)
		
		tt_read_lst += read_lst
		#print(graph_num, sub_graphs_number, tt_graph_num)
	read_lst_lst.sort(key = lambda x:len(x), reverse = True)
	#print('read_lst_lst 1', len(read_lst_lst))
	read_lst_lst = refine_read_lst(read_lst_lst)
	#print('read_lst_lst', len(read_lst_lst))
	o = open('{}/ref_seqs_1st.txt'.format(outdir), 'w')
	o.writelines('')
	o.close()
	read_lst_lst_len = len(read_lst_lst)
	tt_read_seq_dic = extract_reads_from_files(tt_read_lst,  read_chunk_size)
	ref_read_number = 0

	for lst_n in range(int((len(read_lst_lst) - 1) / ((1000 * threads))) + 1):
		read_seq_lst_lst = [[(read_num, tt_read_seq_dic[read_num]) for read_num in read_lst] for read_lst in read_lst_lst[lst_n * 1000 * threads : min(read_lst_lst_len, (lst_n + 1) * 1000 * threads)]]
		#print('len(read_seq_lst_lst)', len(read_seq_lst_lst))
		p = mp.Pool(processes = threads)
		best_seqs_lst_lst = p.map(find_best_seq, read_seq_lst_lst)
		p.close()
		p.join()
		for best_seq_name, best_seq, read_lst in best_seqs_lst_lst:
			if len(read_lst) < 2:
				continue
			ref_seqs_cluster.append(''.join(('{}\t{}'.format(best_seq_name, best_seq), '\t', ''.join([ '{}|{}\t'.format(read_name, score) for read_name, score in read_lst]), '\n')))
			ref_read_number += len(read_lst)
		o = open('{}/ref_seqs_1st.txt'.format(outdir), 'a')		
		o.writelines(ref_seqs_cluster)
		o.close()
		ref_seqs_cluster = []
	print("\033[31m[Find best read finished]\033[0m", 'total {} reads support {} reads, time used: {}  s'.format(ref_read_number, len(sub_graphs), time() - start_time))
	#print('reads supported refs {}'.format(ref_read_number))		
