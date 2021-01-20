CC = g++ 

OMP = 1
FASTQ = 0

DEFINES = -DFASTQ=$(FASTQ) -DOMP=$(OMP)

OMP_LIB = -fopenmp

CPPFLAGS = -Wall -ansi -pedantic -g -O2 -std=c++11 $(DEFINES) $(OMP_LIB)

unBCR_QS_obs = unBCR_QS.o BCRdecode.o 
unBCR_QS: $(unBCR_QS_obs)
	$(CC) $(OMP_LIB) -o unBCR_QS $(unBCR_QS_obs)

clean:
	rm -f core *.o *~ unBCR_QS

depend:
	$(CC) $(OMP_LIB) -MM *.cpp *.c > dependencies.mk

include dependencies.mk
