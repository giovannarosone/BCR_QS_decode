# BCR_QS_decode
This library decodes a collection of sequences from its eBWT computed by BCR. 

Nucleotide sequences datasets often come with per-base quality scores that encode sequencing probability error.

BCR_QS_decode offers two decoding options:

a) **eBWT --> FASTA**

b) **eBWT, QS --> FASTQ**,  
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
./unBCR_QS input output maxLengthRead numthreads
```
where:
- input is  input is the BWT filename with no extension .ebwt (and .ebwt.qs for the QS string)
- output is the output filename
- maxLengthRead is the maximum read length
- numthreads is the maximum number of threads possibly used in parallel regions

### Examples
Decoding with option (a)
```sh
./unBCR_QS test/7seqsVar output_7seqs 13 1
```
we get the file output_7seqs.fasta, which is the string collection in test/7seqsVar.fa.

Decoding with option (b)
```sh
./unBCR_QS test/3seqs output_3seqs 14 1
```
we get the file output_3seqs.fastq, which is the string collection with their quality scores in test/3seqs.fq.
