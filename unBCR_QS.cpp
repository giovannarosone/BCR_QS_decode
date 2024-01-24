//
//  unBCR_QS
//
//  Created by Giovanna on 03/01/2021.
//

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

#include <iostream>
#include <assert.h>
#include <string.h>     // std::string, std::to_string
#include <sstream>
#include <stdio.h>
#include <math.h>

using std::cout;
using std::endl;

#include "BCRdecode.hpp"
#include "Parameters.h"

//#define SIZEBUFFER 1024

using namespace std;


int main(int argc, char *argv[])
{

    #if OMP
        if( argc != 5 ) {
            std::cerr << "usage: " << argv[0] << " input output maxLengthRead numthreads" << std::endl;
    #else
        if( argc != 4 ) {
            std::cerr << "usage: " << argv[0] << " input output maxLengthRead" << std::endl;
    #endif
            std::cerr << "where:" << std::endl;
            std::cerr << "  input is the BWT filename without the extension .ebwt (and .ebwt.qs for the QS string)" << std::endl;
            std::cerr << "  output is the output file " << std::endl;
            std::cerr << "  maxLengthRead is the maximum read length (required)" << std::endl;
            exit(1);
    }

    std::cout << "unBCR_QS: " << argv[0] << std::endl;
    std::cout << "unBCR_QS: The input is " << argv[1] << std::endl;
    std::cout << "unBCR_QS: The output is " << argv[2] << std::endl;
	std::cout << "unBCR_QS: maxLengthRead is " << argv[4] << std::endl;
    if (MODE == 1)
		std::cout << "unBCR_QS: MODE is 1 --> unBCR " << std::endl;
	else if (MODE == 2)
        std::cout << "unBCR_QS: MODE is 2 --> unBCR by using existing partial files" << std::endl;
	else
		std::cout << "unBCR_QS: MODE is 3 --> unBCR by using .table file" << std::endl;
	
    string fileInput=argv[1];
    string fileOutDecode=argv[2];
     
    int num_threads = 1;
    int maxLengthRead = atoi(argv[4]);
    #if OMP
            std::cout << "unBCR_QS: numthreads is " << argv[5] << std::endl;
            num_threads = atoi(argv[5]);
    #endif
    
    //dataTypelenSeq lengthRead=0;    //Length of each text
    //dataTypeNChar lengthTot_plus_eof=0;   //Total length of all texts without $-symbols

    //dataTypeNSeq nText=0;   //number total of texts in filename1
    // dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters

    //dataTypedimAlpha sizeAlpha =0;
    
    
    BCRdecode *BCRdec;
    
    BCRdec = new BCRdecode(fileInput, fileOutDecode, MODE, maxLengthRead, num_threads);
    
    std::cerr << "\nThe rebuilt file is ready! \n";
     
    delete BCRdec;
            
    std::cerr << "The End!\n";

    return 1;
}


