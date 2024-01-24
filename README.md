# BCR_QS_decode
This library decodes a collection of sequences from its eBWT computed by BCR. 

Nucleotide sequences datasets often come with per-base quality scores that encode sequencing probability error.

BCR_QS_decode offers two decoding options:

a) eBWT --> fasta

b) eBWT, QS --> fastq
where QS is the string of quality scores permuted according to the eBWT symbols.

### Install

```sh
git clone https://github.com/giovannarosone/BCR_QS_decode
cd BCR_QS_decode
```

### Compile
One can compile choosing between options (a) and (b) above. 

Decoding option (a) is set by default (FASTQ=0):

```sh
make
```

While for decoding option (b), compile with FASTQ=1:

```sh
make FASTQ=1
```
### Usage

```sh
./unBCR_QS input output mode maxLengthRead numthreads
```
where:
- input is  input is the BWT filename with no extension .ebwt (and .ebwt.qs for the QS string)
- output is the output filename
- maxLengthRead is the maximum read length
- numthreads is the maximum number of threads possibly used in parallel regions
