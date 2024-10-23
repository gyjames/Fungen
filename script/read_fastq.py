#!/home/conda/bin/python

import gzip


def read_fastq_gz(fp):#this code read *.fastq.gz file
	name, seq, qual = '', '', ''
	seq_bool, qual_bool =  False,False
	for line in fp: 
		l = line.decode()
		#print(l[:3], seq_bool, qual_bool)
		if l[0] == '@' and len(seq) <= len(qual): #fastq qual long enough
			if seq != '' and len(seq) == len(qual):#good quality to produce
				yield name, (seq.upper(), qual)
			name = l[1: -1]#remove \n, then start a new one
			seq, qual = '', ''
			seq_bool, qual_bool = True, False
			continue
		if l[0] == '+' and not qual_bool:
			seq_bool, qual_bool = False, True
			continue
		if seq_bool and not qual_bool:
			seq += l.strip()
			continue
		if not seq_bool and qual_bool:
			qual += l.strip()
			continue
	yield name, (seq.upper(), qual)#yield the last read
	yield '####END',('####END', '####END')




def read_fastq(fp): # this is a generator function
	last = None # this is a buffer keeping the last unprocessed line
	while True: # mimic closure; is it a bad idea?
		if not last: # the first record or a record following a fastq
			for l in fp: # search for the start of the next record
				if l[0] in '>@': # fasta/q header line
					last = l[:-1] # save this line
					break
		if not last: break
		name, seqs, last = last[1:].replace(" ", "_"), [], None
		for l in fp: # read the sequence
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1].strip())
		if not last or last[0] != '+': # this is a fasta record
			yield name, (''.join(seqs).upper(), None) # yield a fasta record
			if not last: break
		else: # this is a fastq record
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fp: # read the quality
				seqs.append(l[:-1].strip())
				leng += len(l) - 1
				if leng >= len(seq): # have read enough quality
					last = None
					yield name, (seq.upper(), ''.join(seqs)); # yield a fastq record
					break
			if last: # reach EOF before reading enough quality
				yield name, (seq.upper(), None) # yield a fasta record instead
				break
	yield '####END',('####END', '####END')

def read_fq(file_name):
	if '.gz' in file_name:
		fqs = read_fastq_gz(gzip.open(file_name))
	else:
		fqs = read_fastq(open(file_name))
	return fqs
