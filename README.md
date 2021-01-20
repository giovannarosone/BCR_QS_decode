# BCR_QS_decode
a) eBWT --> FASTA, text

b) eBWT, QS --> FASTQ, text

### Compile
One can compile choosing between options (a) and (b). 

Decoding option (a) (FASTQ=0) is set by default:

```sh
make
```
While for decoding option (b) set FASTQ=1:

```sh
make FASTQ=1
```
### Usage

```sh
./unBCR_QS input output mode maxLengthRead numthreads
```
where:
- input is the filename without extension
- output is the output filename
- mode = 1 --> decode .ebwt just using the maximum read length (maxLengthRead)
- mode = 2 --> decode .ebwt by using existing partial ebwt files, in addition to .info and .table files
- mode = 3 --> decode .ebwt by using .info and .table files
- maxLengthRead is the maximum read length (mandatory for running mode=1)
- numthreads is the maximum number of threads possibly used in parallel regions
