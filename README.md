Fungen is a reference-free tool to construct accurate transcripts using long-read metatranscriptomic data through read clustering and error correction.

Package accessment:

1)	cd path
2)	git clone https://github.com/gyjames/Fungen
3)	cd FunGen

Package install

1)	conda create -n fungen python=3.7
2)	conda activate python=3.7
3)	pip install -r ./requirements.txt

Help infomation

 	python ./FunGen -h

Parameters


	  -h, --help  show this help message and exit 
   
	  -i F   Define the input metatranscriptome fastq file. (default: False)
   
	  -t T   Threads used for clustering. (default: 10)
   
	  -k K   Minimizer size (default: 11)
   
	  -w W   Window size (default: 20)
   
	  -o O   Output folder (default: Recent_dictionary)
   
	  -c C   Read number in each chunk (default: 100000)
   
	  -l ML  Minimum length of compressed read sequence (default: 30)
   
	  -continue U Continue the work from last step using -continue 1 (default: False)
   


Outputs:

          cluster_results.txt:          The final clustering results

          cluster_represents.fasta:     Reference sequences for gene clusters

          corrected_seqs.fastq:         Corrected reads

