# BCR_QS_decode
a) eBWT --> FASTA, text

b) eBWT, QS --> FASTQ, text

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
- input is the filename without extension
- output is the output filename
- mode can be set equal to
  * 1 --> for decoding eBWT by using the maximum read length only (maxLengthRead)
  * 2 --> for decoding eBWT by using existing partial ebwt files, in addition to .info and .table files
  * 3 --> for decoding eBWT by using .info and .table files
- maxLengthRead is the maximum read length (mandatory for running mode=1)
- numthreads is the maximum number of threads possibly used in parallel regions
