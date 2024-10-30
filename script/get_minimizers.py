#!/usr/bin/python3
import os,sys
from collections import deque, Counter
import itertools
import multiprocessing as mp
from time import time
from script.read_fastq import read_fq
import random


def rev_string(string):
	return ''.join([string[len(string)-n-1] for n in range(0,len(string))])

def rev_com_seq(sequence):
	nuc_rev={'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a','N':'N','X':'X','n':'n','Y':'R','R':'Y','K':'M','M':'K','S':'S','W':'W','B':'V','V':'B','H':'D','D':'H','y':'r','r':'y','k':'m','m':'k','s':'s','w':'w','b':'v','v':'b','h':'d','d':'h', '0': '3', '3': '0', '2': '1', '1': '2', '4': '4'} 
	return ''.join([nuc_rev[sequence[len(sequence)-n-1]] for n in range(0,len(sequence))])
 
def mul_reduce(l):
	lr=1
	for lx in l:
		lr*=lx
	return lr


def poisson_mean(qual_str):
	D = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)}
	qual = [q for q in qual_str]
	if len(qual) == 0:
		return 0
#       return sum([qual.count(char_)*(1-D[char_]) for i,char_ in enumerate(D) ])/len(qual)
	#print(sum([qual.count(char_)*(1-D[char_]) for i,char_ in enumerate(D) ])/len(qual),sum(((1-D[char_]) * n for char_,n in Counter(qual).most_common()))/len(qual))
	return sum(((1-D[char_]) * n for char_,n in Counter(qual).most_common()))/len(qual)


def compress_sequence_bak(sequence,qual):
	compressed_sequence = ''.join(chr_ for chr_,_ in itertools.groupby(sequence))
	compressed_num = [len([z for z in ob]) for chr_,ob in itertools.groupby(sequence)]
	de_qual,compressed_qual = deque(qual),[]
	for k in compressed_num:
		compressed_qual.append(min([de_qual[kn] for kn in range(k)]))
		t = [de_qual.popleft() for i in range(0,k)]
	compressed_num_str = ''.join([str(k)+':' for k in compressed_num])
	compressed_qual_str = ''.join([str(q) for q in compressed_qual])
	return compressed_sequence,compressed_qual_str,compressed_num_str

def compress_sequence_qual(sequence,qual):
	nc2num = {'A':'0', 'a':'0', 'C':'1', 'c':'1', 'G':'2', 'g':'2', 'T':'3', 't':'3', 'N': '4', 'n': '4'}
	compressed_sequence = ''.join([nc2num[chr_] for chr_,_ in itertools.groupby(sequence) if chr_ in nc2num])
	compressed_num = [len([z for z in ob]) for chr_,ob in itertools.groupby(sequence)]
	de_qual,compressed_qual = deque(qual),[]
	for k in compressed_num:
		compressed_qual.append(min([de_qual[kn] for kn in range(k)]))
		t = [de_qual.popleft() for i in range(0,k)]
	compressed_num_str = ''.join([str(k)+':' for k in compressed_num])
	compressed_qual_str = ''.join([str(q) for q in compressed_qual])
	return compressed_sequence,compressed_qual_str,compressed_num_str

def compress_sequence(sequence):
	nc2num = {'A':'0', 'a':'0', 'C':'1', 'c':'1', 'G':'2', 'g':'2', 'T':'3', 't':'3', 'N': '4', 'n': '4'}
	compressed_sequence = ''.join([nc2num[chr_] for chr_,_ in itertools.groupby(sequence) if chr_ in nc2num])
	compressed_num = [len([z for z in ob]) for chr_,ob in itertools.groupby(sequence)]
	#compressed_num_str = ''.join([str(k)+':' for k in compressed_num])
	return compressed_sequence, compressed_num


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
	return minimizers_hash	

def get_minimizers_seq_qual(sequence, qual, k, w):
	if len(sequence) <= k:
		return []
	global D
	D = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)}
	minimizers = deque([(sequence[i: i+k], i) for i in range(w - k + 1)])
	minimizer = min(minimizers)
	minimizer_seq,minimizer_lo = minimizer
	minimizer_qual = poisson_mean(qual[minimizer_lo: minimizer_lo+k])
	minimizers_hash = [(minimizer_seq,minimizer_qual, minimizer_lo)]
	for i in range(w - k + 1,len(sequence) - k):
		discard_kmer, new_kmer = minimizers.popleft(),(sequence[i: i+k], i)
		minimizers.append(new_kmer)
		if discard_kmer == minimizer:
			minimizer = min(minimizers)
			minimizer_seq, minimizer_lo = minimizer
			minimizer_qual = minimizer_poisson_mean(qual[minimizer_lo: minimizer_lo+k])
			new_minimizer = (minimizer_seq,minimizer_qual, minimizer_lo)
			minimizers_hash.append(new_minimizer)
		elif new_kmer < minimizer:
			minimizer = new_kmer
			minimizer_seq, minimizer_lo = minimizer
			minimizer_qual = minimizer_poisson_mean(qual[minimizer_lo: minimizer_lo+k])
			new_minimizer = (minimizer_seq, minimizer_qual, minimizer_lo)
			minimizers_hash.append(new_minimizer)
	return minimizers_hash


def generate_minimizers(seq_info):
	seq_num, seq_name_raw, seq, k, w = seq_info
	seq_rev = rev_com_seq(seq)
	compressed_seq, compress_num = compress_sequence(seq)
	compress_num_str = ''.join([str(k)+':' for k in compress_num])
	compressed_seq_rev = rev_com_seq(compressed_seq)
	compress_numv = reversed(compress_num)
	compress_num_strv = ''.join([str(k)+':' for k in compress_numv])
	query_minimizers_raw = [minimizer for minimizer, i in get_minimizers_seq(compressed_seq, k, w)]
	query_minimizers = [int(minimizer) for minimizer in query_minimizers_raw]
	query_minimizers_rev = [int(minimizer) for minimizer, i in get_minimizers_seq(compressed_seq_rev,k, w)]
	if len(query_minimizers) < 10 or len(query_minimizers_rev) < 10:
		query_minimizers, query_minimizers_rev = [], []
	elif len(list(set(query_minimizers))) / len(query_minimizers) < 0.5 or len(list(set(query_minimizers_rev))) / len(query_minimizers_rev) < 0.5:
		query_minimizers, query_minimizers_rev = [], []
	if len(list(set(query_minimizers))) < len(query_minimizers) / 2: #remove simple repeats
		query_minimizers, query_minimizers_rev = [], []
	if query_minimizers != []:
		compressed_num = [len([z for z in ob]) for chr_,ob in itertools.groupby(query_minimizers)]
		compressed_num.sort(reverse = True)
		if compressed_num[0] > 20:
			query_minimizers, query_minimizers_rev = [], []
	#print('compressed_seq, compress_num_str, compressed_seq_rev, compress_num_strv', compressed_seq, compress_num_str, compressed_seq_rev, compress_num_strv)
	#return (seq_num, seq_name_raw, query_minimizers, query_minimizers_rev, compressed_seq, compress_num_str, compressed_seq_rev, compress_num_strv)
	return (seq_num, seq_name_raw, query_minimizers, query_minimizers_rev, compressed_seq, seq, compressed_seq_rev, seq_rev)

def sort_minimizer_name(outdir):
	minimizer_name_number_data = []
	with open('{}/name_num_record.txt'.format(outdir)) as name_num_record:
		for line in name_num_record:
			seq_num, seq_name, minimizer_num = line.split('\t')
			minimizer_name_number_data.append((int(seq_num), int(minimizer_num.strip())))
	minimizer_name_number_data.sort(key = lambda x:x[1], reverse = True)
	name_num_record_sort_file = open('{}/name_num_record_sort.txt'.format(outdir), 'w')
	for seq_num, minimizer_num in minimizer_name_number_data:
		name_num_record_sort_file.writelines('{}\t{}\n'.format(seq_num, minimizer_num))
	name_num_record_sort_file.close()
	
	sorted_read_number =  len(minimizer_name_number_data)
	del minimizer_name_number_data
	return sorted_read_number

def get_minimizers_from_file(fqs_name, k, w, threads, chunk_read_num, outdir):
	os.mkdir('{}/minimizer_record'.format(outdir))
	os.mkdir('{}/seq_record'.format(outdir))
	name_record_data = []
	minimizer_bit_dic = {}
	minimizer_all_lst = [[]]
	minimizer_num_lst = []
	minimizer_name = []
	seq_record_data = []
	minimizer_num = 0
	name_record_file = open('{}/name_num_record.txt'.format(outdir), 'w')
	name_record_file.writelines('')
	name_record_file.close()
	minimizer_set = set([])
	minimizer_counter = Counter()
	for t,(name,(seq, qual)) in enumerate(read_fq(fqs_name)):
		seq_name = t + 1
		seq_namev = -t - 1
		minimizer_all_lst.append([])
		if name != '####END':
			minimizer_name.append((seq_name, name, seq, k, w))
		if (t + 1) % chunk_read_num == 0 or (name == '####END' and seq == '####END'):
			read_time = time()
			record_data = []
			seq_record_data = []
			chunk_name = int(t / chunk_read_num)
			p = mp.Pool(processes = threads)
			minimizer_lst = p.map(generate_minimizers, minimizer_name)
			p.close()
			p.join()
			minimizer_name = []
			minimizer_record = []
			#for seq_name, seq_name_raw, query_minimizers, query_minimizers_rev, compressed_seq, compress_num_str, compressed_seq_rev, compress_num_strv in minimizer_lst:
			for seq_name, seq_name_raw, query_minimizers, query_minimizers_rev, compressed_seq, seq, compressed_seq_rev, seq_rev in minimizer_lst:
				if len(query_minimizers) <= 4 or len(query_minimizers) > 2000:
					minimizer_str = ''.join((str(seq_name), '\t', '-', '\n'))
					record_data.append(minimizer_str)
					seq_str = ''.join((str(seq_name), '\t', '-', '\n'))
					seq_record_data.append(seq_str)
					continue
				minimizer_set.update(query_minimizers)
				minimizer_set.update(query_minimizers_rev)
				if chunk_name % 10 == 0:
					minimizer_counter.update(query_minimizers)
					minimizer_counter.update(query_minimizers_rev)
				name_record_data.append('{}\t{}\t{}\n'.format(seq_name, seq_name_raw, len(query_minimizers)))
				minimizer_str = ''.join((str(seq_name), '\t', ''.join(([str(w) + ' ' for w in query_minimizers + ['|'] + query_minimizers_rev])), '\n'))
				record_data.append(minimizer_str)
				#seq_str = ''.join((str(seq_name), '\t', compressed_seq, '\t', compress_num_str, compressed_seq_rev, '\t', compress_num_strv, '\n'))
				seq_str = ''.join((str(seq_name), '\t', compressed_seq, '\t', seq, '\t', compressed_seq_rev, '\t', seq_rev, '\n'))
				seq_record_data.append(seq_str)
			record_data_file = open('{}/minimizer_record/{}'.format(outdir, chunk_name), 'w')
			record_data_file.writelines(record_data)
			record_data_file.close()
			seq_data_file = open('{}/seq_record/{}'.format(outdir, chunk_name), 'w')
			seq_data_file.writelines(seq_record_data)
			seq_data_file.close()
			minimizer_num = len(list(minimizer_set))
			print("\033[31m[Call minimizer]\033[0m",'processed_read: {};'.format(t + 1), 'identified_minimizers: {};'.format(minimizer_num), 'time_used: %.3f s'%( time() - read_time))
			record_data = []
			del record_data
			name_record_file = open('{}/name_num_record.txt'.format(outdir), 'a')
			name_record_file.writelines(name_record_data)
			name_record_file.close()
			name_record_data = []
			seq_record_data = []
			del seq_record_data
	minimizer_counter.update(list(minimizer_set))
	minimizer_lst = [''.join((str(minimizer), '\t', str(d))) for minimizer, d in list(minimizer_counter.most_common())]
	#minimizer_lst = list(minimizer_set)
	minimizer_lst_record_data = ''.join([''.join((str(minimizer_lst[i]), '\t', str(i), '\n')) for i in range(len(minimizer_lst))])
	minimizer_lst_record_file = open('{}/minimizer_lst_array.div'.format(outdir), 'w')
	minimizer_lst_record_file.writelines(minimizer_lst_record_data)
	minimizer_lst_record_file.close()


def decord_minimizers(record_data):
	num_data, minimizer_data = record_data.split('\t')[:2]
	read_num = int(num_data)
	read_numv = -int(num_data)
	query_minimizer_data = minimizer_data.split('|')[0].strip()
	query_minimizer_data_rev = minimizer_data.split('|')[1].strip()
	query_minimizers = [int(w) for w in query_minimizer_data.split()]
	query_minimizers_rev = [int(w) for w in query_minimizer_data_rev.split()]
	return read_num, read_numv, query_minimizers, query_minimizers_rev



def extract_sorted_minimizers_from_file(reads_in_chunk, t, chunk_read_num, outdir):
	#print(len(reads_in_chunk))
	extract_time = time()
	chunk_name = int((t - 1) / chunk_read_num)
	reads_file_dic = {}
	chunk_data = []
	seq_chunk_data = []
	for seq_num, order_t in reads_in_chunk:
		file_name, record_num = int((seq_num -1) / chunk_read_num), (seq_num - 1) % chunk_read_num + 1
		if file_name not in reads_file_dic:
			reads_file_dic[file_name] = []
		reads_file_dic[file_name].append((record_num, order_t))
	tt = 0
	reads_file_lst = []
	for i,file_name in enumerate(reads_file_dic):
		reads_file_lst.append((file_name, reads_file_dic[file_name]))
	random.shuffle(reads_file_lst)

	max_minimizer_num = 0
	min_minimizer_num = 10 ** 10

	for file_name,read_record_lst in reads_file_lst:
		#minimizer_data = open('minimizer_record/{}'.format(file_name)).readlines()
		read_record_lst.sort(key = lambda x:x[0])
		#print(read_record_lst[:10])
		read_record_deque = deque(read_record_lst)
		tt = 0
		if len(read_record_deque) != 0:
			fst_read_num, fst_read_order = read_record_deque.popleft()
		else:
			break
		line_num = 0
		with open('{}/minimizer_record/{}'.format(outdir, file_name)) as minimizer_record_file:
			for record_line in minimizer_record_file:
				line_num +=1
				if line_num == fst_read_num:
					tt += 1
					read_num, read_numv, query_minimizers, query_minimizers_rev = decord_minimizers(record_line)
					if len(query_minimizers) > max_minimizer_num:
						max_minimizer_num = len(query_minimizers)
					if len(query_minimizers) < min_minimizer_num:
						min_minimizer_num = len(query_minimizers)
					query_minimizers_bit = [MINIMIZER_BIT_DIC[minimizer] for minimizer in query_minimizers]
					query_minimizers_bit_rev = [MINIMIZER_BIT_DIC[minimizer] for minimizer in query_minimizers_rev]
					record_line_new = ''.join((str(read_num), '\t', ''.join(([str(w) + ' ' for w in query_minimizers_bit + ['|'] + query_minimizers_bit_rev])), '\n'))
					chunk_data.append((record_line_new, fst_read_order))
					if len(read_record_deque) == 0:
						break
					else:
						fst_read_num, fst_read_order = read_record_deque.popleft()


		read_record_deque = deque(read_record_lst)
		if len(read_record_deque) != 0:
			fst_read_num, fst_read_order = read_record_deque.popleft()
		else:
			break
		line_num = 0
		with open('{}/seq_record/{}'.format(outdir, file_name)) as seq_record_file:
			for seq_record_line in seq_record_file:
				line_num +=1
				if line_num == fst_read_num:
					compressed_seq = ''.join((seq_record_line.split('\t')[0], '\t', seq_record_line.split('\t')[1], '\t', seq_record_line.split('\t')[3]))
					seq_chunk_data.append((''.join((compressed_seq, '\n')), fst_read_order))
					if len(read_record_deque) == 0:
						break
					else:
						fst_read_num, fst_read_order = read_record_deque.popleft()
	
	chunk_data.sort(key = lambda x:x[1])
	seq_chunk_data.sort(key = lambda x:x[1])
	chunk_output_data = [read_minimizer_line for read_minimizer_line, order_t in chunk_data]
	seq_chunk_output_data = [seq_line for seq_line, order_t in seq_chunk_data]
	output_file = open('{}/minimizer_record_sort/{}'.format(outdir, chunk_name), 'w')
	output_file.writelines(chunk_output_data)
	output_file.close()
	output_len_file = open('{}/minimizer_record_len/{}'.format(outdir, chunk_name), 'w')
	output_len_file.writelines('{}\t{}\n'.format(max_minimizer_num, min_minimizer_num))
	output_len_file.close()
	output_seq_file = open('{}/seq_record_sort/{}'.format(outdir, chunk_name), 'w')
	output_seq_file.writelines(seq_chunk_output_data)
	output_seq_file.close()

	print("\033[31m[Call minimizer]\033[0m",'sorted_read_chunk: {}-{};'.format(t - chunk_read_num, t), 'time_used: %.3f s'%( time() - extract_time))

def initializer_extract_sorted_minimizers_from_file(minimizer_bit_dic):
	global MINIMIZER_BIT_DIC
	MINIMIZER_BIT_DIC = minimizer_bit_dic


def sort_minimizers_from_file(chunk_read_num, sorted_read_number, threads, outdir):
	sort_time = time()
	os.mkdir('{}/minimizer_record_sort'.format(outdir))
	os.mkdir('{}/minimizer_record_len'.format(outdir))
	os.mkdir('{}/seq_record_sort'.format(outdir))
	reads_in_chunk = []
	minimizer_bit_dic = {}
	with open('{}/minimizer_lst_array.div'.format(outdir)) as minimizer_lst_file:
		for line in minimizer_lst_file:
			minimizer, d, minimizer_bit = line.split()
			minimizer_bit_dic[int(minimizer)] = int(minimizer_bit)
	p = mp.Pool(processes = min(threads, 5), initializer = initializer_extract_sorted_minimizers_from_file(minimizer_bit_dic))
	t = 0
	with open('{}/name_num_record_sort.txt'.format(outdir)) as name_num_record_sort_file:
		for line in name_num_record_sort_file:
			seq_num = int(line.split('\t')[0])
			t += 1
			reads_in_chunk.append((seq_num, t))
			if t % chunk_read_num == 0 or t == sorted_read_number:
				#extract_sorted_minimizers_from_file(reads_in_chunk, t, chunk_read_num)
				p.apply_async(extract_sorted_minimizers_from_file, (reads_in_chunk, t, chunk_read_num, outdir, ))
				reads_in_chunk = []
	p.close()
	p.join()
			
		

def get_minimizers(args, outdir):
	print("\033[31m**************\033[0m", 'Extract minimizers from files', "\033[31m**************\033[0m")

	try: os.mkdir('{}'.format(outdir))
	except:
		os.system("rm -rf {}".format(outdir))
		os.mkdir('{}'.format(outdir))
		
	

	get_minimizer_start_time = time()
	k = args.k
	w = args.w
	threads = args.t
	chunk_read_num = args.c
	fqs_name = args.f
	get_minimizers_from_file(fqs_name, k, w, threads, chunk_read_num, outdir)
	read_time = time()
	sorted_read_number = sort_minimizer_name(outdir)
	print("\033[31m[Call minimizer]\033[0m",'sorted reads: {};'.format(sorted_read_number), 'time_used: %.3f s'%( time() - read_time))
	sort_minimizers_from_file(chunk_read_num, sorted_read_number, threads, outdir)
	print("\033[31m[Call minimizer finished]\033[0m",'Total sorted reads: {};'.format(sorted_read_number), 'time_used: %.3f s'%( time() - get_minimizer_start_time))
	


#if __name__=='__main__':
#	get_minimizers()
#	print(compress_sequence('aaaaaannaaaa'))
#	print('AGACCATGACTATA',get_minimizers('AGACCATGACTA',13,15))
