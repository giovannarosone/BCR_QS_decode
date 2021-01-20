//
//  unBCRforPosition.hpp
//  unBCRforPosition
//
//  Created by Giovanna on 27/12/20.
//

#ifndef BCRdecode_hpp
#define BCRdecode_hpp

/**
 ** BCR is part of:
 ** BEETL: Burrows-Wheeler Extended Tool Library
 ** Documentation in: doc/BEETL.md
 **
 ** Copyright (c) 2011-2014 Illumina, Inc. **
 ** BEETL software package is
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citations:
 **
 ** Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 **
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone and Marinella Sciortino
 ** Lightweight LCP Construction for Next-Generation Sequencing Datasets.
 ** Proceedings of WABI 2012, pp 326-337, 2012
 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone
 ** Lightweight algorithms for constructing and inverting the BWT of string collections.
 ** Theoretical Computer Science 483: 134-148 (2013)
 **
 ** Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
 ** Lightweight LCP construction for very large collections of strings.
 ** Journal of Discrete Algorithms (2016) 37: 17-33
 **
 ** By Giovanna Rosone
 **
 **/


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
//#include <stdlib.h>

#include "Parameters.h" // Defines ulong and uchar.
#include <string.h>

#include <deque>


#define BUFFERSIZE 1024

#define SIZE_ALPHA 256

#define DIMBLOCK  1024     // 1048576

using std::string;
using std::vector;

class BCRdecode
{
public:
    BCRdecode(string, string, int, dataTypelenSeq, int);
    ~BCRdecode();

    string  fileOutBwt;
    string  fileOutCyc;
    string  ext;
    
    struct sortElement {
        dataTypeNSeq seqN;
        dataTypeNChar posN;
    };
    
    dataTypeNChar SIZEBUFFERcycFiles;
    
    dataTypelenSeq lengthRead;    //Length of each text
    dataTypeNChar lengthTot;   //Total length of all texts without $-symbols
    dataTypeNChar lengthTot_plus_eof;   //Total length of all texts with $-symbols

    dataTypeNSeq nText;   //number total of texts in filename1
    dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters
    
    int numthreads;
    
    dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
    dataTypedimAlpha alpha[SIZE_ALPHA]; //Corresponding between the alphabet, the piles and tableOcc
    dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
    dataTypedimAlpha *alphaInverse;  //Corresponding between alpha[i] and the symbol as char
    
    std::vector< std::deque<sortElement> > vectPair;  //Is is used both encoding, decoding, searching.
            //ulong seqN;  //contains a number of a sequence
            //ulong posN;  //contains the position of the last inserted symbol of the sequence seqN[i]
            //uchar pileN; //contains the number of the pile of the last inserted symbol of the sequence seqN[i]
    
    std::vector< std::deque<sortElement> > vectPairTmp;
    std::vector< std::vector< std::deque<sortElement> > > vectPairParallel;
    std::vector< std::vector< std::deque<sortElement> > > vectPairParallelTmp;
    
    std::vector <bool> vectInsTexts;
    dataTypeNSeq textToInsert;
    vector <dataTypeNChar> vectSizeCurrentPile;

    vector< vector< vector<dataTypeNChar> > > vectorOcc;
    vector <dataTypeNChar> numBlocksInPartialBWT;

    int concatenateFastaOrFastq( string );
    
    int recoverInfo(string , int );

    int splitIntoPartial(string, int);
    
    int buildFreq(string);
  
    int decodeBCRmultipleReverse(string);
    #if OMP
        //16-01
        //int RecoverNsymbolsReverseByVector(int t_id, string file1, uchar * newSymb, uchar * newQual);
        int RecoverNsymbolsReverseByVector(int t_id, string file1, uchar * newSymb);
    #else
        //16-01
        //int RecoverNsymbolsReverseByVector(string file1, uchar * newSymb, uchar * newQual);
        int RecoverNsymbolsReverseByVector(string file1, uchar * newSymb);
    #endif
    
    dataTypeNChar rankManySymbolsByVector(FILE & , dataTypeNChar *, dataTypeNChar, uchar *, uchar *, FILE &);
    
    int computeVectorUnbuildBCR();
    
    int convertFromCycFileToFastaOrFastq( string );

    
private:

    int findBlockToRead(dataTypeNChar *, dataTypedimAlpha , dataTypeNChar *, dataTypeNChar *);
    
    FILE * openFilePartialIn(dataTypedimAlpha currentPile);
    dataTypeNChar readOnFilePartial(uchar *buffer, dataTypeNChar toRead, FILE * InFileBWT);
    dataTypeNChar writeOnFilePartial(uchar *buffer, dataTypeNChar numchar, FILE * OutFileBWT);
    int closeFilePartial(FILE * pFile);


    
};




#endif /* BCRdecode_hpp */
