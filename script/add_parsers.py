#!/home/conda/bin/python
import argparse

def add_parsers():
	parsers  =  argparse.ArgumentParser(description = "Reference-free clustering of long-read fungal metatranscriptome reads", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
	#print('Basic arguments:\n','-'*30)	
	parsers.add_argument('-i', dest = "f",type = str,  default = False, help = 'Define the input metatranscriptome fastq file.')
	parsers.add_argument('-t', dest = "t", type = int, default = 10, help = 'Threads used for clustering.')
	parsers.add_argument('-k', type = int, default = 11, help = 'Minimizer size')
	parsers.add_argument('-w', type = int, default = 20, help = 'Window size')
	parsers.add_argument('-o', type = str, default = 'Recent_dictionary', help = 'Output folder')
	parsers.add_argument('-c', type = int, default = 10 ** 5, help = 'Read number in each chunk')
	parsers.add_argument('-l', dest = 'ml', type = int, default = 30, help = 'Minimum length of compressed read sequence')
	parsers.add_argument('-continue', dest = 'u', type = bool, default = False, help = 'Continue the work from last step using -continue 1')
	return parsers

