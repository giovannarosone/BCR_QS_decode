# BCR_QS_decode
eBWT --> FastQ, FASTA, text

usage: ./unBCR_QS input output mode maxLengthRead numthreads
where:
  input is the filename without extension
  output is the output filename
	mode = 1 --> decode .ebwt just using the maximum read length (maxLengthRead)
	mode = 2 --> decode .ebwt by using existing partial ebwt files, in addition to .info and .table files
	mode = 3 --> decode .ebwt by using .info and .table files
  maxLengthRead is the maximum read length (mandatory for running mode=1)
