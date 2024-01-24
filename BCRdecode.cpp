//
//  BCRdecode.cpp
//  BCRdecode
//
//  Created by Giovanna on 04/01/2021.

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
 
#include "BCRdecode.hpp"
//#include "Tools.h"
#include "Parameters.h"

#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <string.h>     // std::string, std::to_string
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>      // std::stringstream

#include <math.h>

#include <omp.h>

#include <sys/resource.h>

using namespace std;


BCRdecode::BCRdecode (string fileInput, string fileOutDecode, int mode, dataTypelenSeq maxLength, int num_threads)
{
    
    fileOutBwt = "bwt_";
    char c_aux[500];
    strcpy (c_aux,"mkdir -p ");
    strcat (c_aux, "cycFilesDecode/");
    assert (system (c_aux) == 0);
    fileOutCyc = "cycFilesDecode/cyc.\0";
    
    ext = ".aux";

    lengthRead=maxLength;
    
    #if OMP
        numthreads = num_threads; //number of threads set by the user
    #endif
    
    if( mode==2 || mode==3 ) //recoverInfo
    {
        nText=0;  //number total of texts in fileInput
        lengthTot_plus_eof=0; //length of the BWT
        sizeAlpha=0;  //number of the different symbols in the input texts
        
        assert ( recoverInfo(fileInput, mode) == 1);
        
        #if OMP
            if((int)nText<numthreads){
                std::cerr << "The number of sequences is " << nText << ", which is too small to use " << numthreads << " threads.\n Number of threads used: " << nText << std::endl;
                numthreads = nText;
            }
        #endif
        
        if (mode == 3) {
            time_t startS,endS;
            double difS;
            
            time (&startS);
            //run splitIntoPartial to compute BWT-partial
            std::cerr << "Start splitIntoPartial " << startS << " seconds\n";
            assert (splitIntoPartial(fileInput,3) == 1);
            time (&endS);
            difS = difftime (endS,startS);
            std::cerr << "End splitIntoPartial " << endS << " seconds\n";
            std::cerr << "splitIntoPartial tooks " << difS << " seconds\n";
        }
        else
            std::cerr << "\nWe skip splitIntoPartial, we have BWT/QS partial files.\n";
    }
    
    if (mode == 1)
    {
        time_t startS,endS;
        double difS;
        //run buildFreq
        assert (buildFreq(fileInput) == 1);
        #if OMP
            if((int)nText<numthreads){
                std::cerr << "The number of sequences is " << nText << ", which is too small to use " << numthreads << " threads.\n Number of threads used: " << nText << std::endl;
                numthreads = nText;
            }
        #endif
        time (&startS);
        //run splitIntoPartial to compute BWT-partial
        std::cerr << "Start splitIntoPartial " << startS << " seconds\n";
        assert (splitIntoPartial(fileInput,1) == 1);
        time (&endS);
        difS = difftime (endS,startS);
        std::cerr << "End splitIntoPartial " << endS << " seconds\n";
        std::cerr << "splitIntoPartial tooks " << difS << " seconds\n";

        #if ((DEBUG == 1) || (verboseDecode == 1))
            std::cerr << "i" << "\t" << "freq" << "\t" << "Code" << "\t" << "ASCII" << "\n";
            dataTypedimAlpha mmm=0;
            for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; i++)
                if (freq[i] > 0) {
                    std::cerr << i << "\t" << freq[i] << "\t" << (unsigned int)alpha[i] << "\t" << (unsigned int)alphaInverse[mmm] << "\n";
                    mmm++;
                }
        
            std::cerr << "TableOcc: "  << "\n";
            for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
                std::cerr << int(g)  << ":\t";
                for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++)
                    std::cerr << tableOcc[g][j]  << "\t";
                std::cerr << "\n";
            }
        #endif
    }
    
    
    std::cerr << "For the computation of the new positon useful for unbuildBCR we use the vector of the occurrences of " << DIMBLOCK << " size" << std::endl;
    
    time_t startC,endC;
    double difC;
    time (&startC);
    std::cerr << "\nStart computeVectorUnbuildBCR " << startC << " seconds\n";
    
    //run computeVectorUnbuildBCR to set vectorOcc
    
    assert ( computeVectorUnbuildBCR() == 1);

    time (&endC);
    difC = difftime (endC,startC);
    std::cerr << "End    computeVectorUnbuildBCR " << endC << " seconds\n";
    std::cerr << "computeVectorUnbuildBCR tooks " << difC << " seconds\n";

#if DEBUG
    for (int x = 0 ; x < sizeAlpha; x++)  {
        for(int z = 0; z < sizeAlpha ; z++)   {
            cout << "vectorOcc[" << x << "][" << z << "]=" << vectorOcc[x][z][0] << "\t";
        }
        cout << endl;
    }
#endif

    time_t startI,endI;
    double difI;
    time (&startI);
    std::cerr << "Inverse BWT by Backward direction." << std::endl;
    std::cerr << "\nStart decodeBCRmultipleReverse " << startI << " seconds\n";

    //run decodeBCRmultipleReverse
    textToInsert = nText;
    assert ( decodeBCRmultipleReverse(fileInput) == 1);
    std::cerr << "The cyc files have been built in "<< "cycFilesDecode" <<" folder. Building the sequences." << std::endl;
    
    time (&endI);
    difI = difftime (endI,startI);
    std::cerr << "End   decodeBCRmultipleReverse " << endI << " seconds\n";
    std::cerr << "decodeBCRmultipleReverse tooks " << difI << " seconds\n";

    string newFilename;
    #if USE_QS == 1
        newFilename = fileOutDecode + ".fastq";
    #else
        newFilename = fileOutDecode + ".fasta";
    #endif


    time_t startT,endT;
    time (&startT);
    
    std::cerr << "\nStart convertFromCycFileToFastaOrFastq " << startT << " seconds\n";
    
    //run convertFromCycFileToFastaOrFastq
    assert (convertFromCycFileToFastaOrFastq( newFilename) == 1);
    
    time (&endT);
    std::cerr << "End   convertFromCycFileToFastaOrFastq " << endT << " seconds\n";
    std::cerr << "convertFromCycFileToFastaOrFastq tooks " << difftime (endT,startT) << " seconds\n";
    
    #if OMP
    //concatenate multiple FastaOrFastq files: newFilename_1....newFilename_(numthread-1) in newFilename
    if(numthreads>1){
        time (&startT);
        std::cerr << "\nStart concatenateFastaOrFastq " << startT << " seconds\n";
        assert(concatenateFastaOrFastq(newFilename)==1);
        time (&endT);
        std::cerr << "End   concatenateFastaOrFastq " << endT << " seconds\n";
        std::cerr << "concatenateFastaOrFastq tooks " << difftime (endT,startT) << " seconds\n";
    }
    #endif
    
    //Free the memory
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        delete [] tableOcc[j];
        tableOcc[j] = NULL;
    }
    delete [] tableOcc;

    delete[] alphaInverse;
    
    /////////////
    char *filenameIn = new char[120];
    
    #if (deleteCycFiles == 1)
        for (dataTypelenSeq g = 0 ; g < lengthRead; g++) {
            sprintf (filenameIn,"%s%d%s",fileOutDecode,g,".txt");
            if (remove(filenameIn)!=0)
                std::cerr << "unbuildBCR: Error deleting " << filenameIn << " file" << std::endl;
        }
        #if ((deletePartialQS == 1) && (USE_QS==1))
            for (dataTypelenSeq g = 0 ; g < lengthRead; g++) {
                sprintf (filenameIn,"%sqs.%d%s",fileOutDecode,g,".txt");
                if (remove(filenameIn)!=0)
                    std::cerr << "unbuildBCR: Error deleting " << filenameIn << " file" << std::endl;
            }
        #endif
    #endif
    
    #if (deletePartialBWT == 1) && (KEEP_eBWT_IN_EXT_MEMORY == 1) || ((deletePartialQS == 1) && (USE_QS==1))
        //|| ((deletePartialLCP == 1) && (BUILD_LCP == 1))
        // ||  (deletePartialGSA == 1 && (BUILD_DA == 1) )
        // || ( deletePartialGSA == 1 && (BUILD_SA == 1) )
        char *filename = new char[100];
        for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
            #if (deletePartialBWT == 1) && (KEEP_eBWT_IN_EXT_MEMORY == 1)
                sprintf (filename, "bwt_%d", g);
                sprintf (filenameIn,"%s%s",filename,ext);
                if (remove(filenameIn)!=0)
                    std::cerr << "unbuildBCR: Error deleting bwt aux files" << std::endl;
            #endif

            #if ((deletePartialQS == 1) && (USE_QS==1))
                sprintf (filename, "qs_%d", g);
                sprintf (filenameIn,"%s%s",filename,ext);
                if (remove(filenameIn)!=0)
                        std::cerr << "unbuildBCR: Error deleting QS aux file" << std::endl;
            #endif

            /*
            #if ((deletePartialLCP == 1) && (BUILD_LCP == 1))
                sprintf (filename, "lcp_%d", g);
                sprintf (filenameIn,"%s%s",filename,ext);
                if (remove(filenameIn)!=0)
                    std::cerr << "unbuildBCR: Error deleting lcp aux files" << std::endl;
            #endif

            #if (deletePartialGSA == 1 && (BUILD_DA == 1) )
                sprintf (filename, "da_%d", g);
                sprintf (filenameIn,"%s%s",filename,ext);
                if (remove(filenameIn)!=0)
                    std::cerr << "unbuildBCR: Error deleting da aux files" << std::endl;
            #endif

            #if deletePartialGSA == 1 && (BUILD_SA == 1)
                sprintf (filename, "sa_%d", g);
                sprintf (filenameIn,"%s%s",filename,ext);
                if (remove(filenameIn)!=0)
                    std::cerr << "unbuildBCR: Error deleting sa aux files" << std::endl;
            #endif
             */
        }
        delete [] filename;
    #endif
    delete [] filenameIn;
   

    
}

//14-01-2021
int BCRdecode::concatenateFastaOrFastq (string Filename) {
    
    string toCopyFile;
    int ret;
    
    //Open newFilename to append
    std::ofstream firstFile;
    firstFile.open(Filename.c_str(), std::ios_base::app);
    
    for ( int i = 1; i<numthreads ; i++){
        //Open newFilename_i
        toCopyFile = Filename + "_" + to_string(i);
        std::ifstream nextFile(toCopyFile.c_str()); // Open for reading
        
        string line;
        
        while(getline(nextFile,line))//Append content to newFilename
            firstFile << line << "\n";
        
        nextFile.close();
        
        //Remove it
        ret = remove(toCopyFile.c_str());
        if(ret == 0)
            std::cout << "File " << toCopyFile << " deleted successfully\n";
        else
            std::cout << "Error: unable to delete file " << toCopyFile << std::endl;;
    }//end-for
    firstFile.close();
    
    return 1;
}


int BCRdecode::recoverInfo(string filenameInfo, int mode) {           //
    

    string fnAuxBCR = filenameInfo + ".info\0";
    FILE* InfoFileBCR = fopen(fnAuxBCR.c_str(), "rb");
    if (InfoFileBCR==NULL) {
        std::cerr << "(lengthBWT+NSequences+sizeAlpha+maxLenSequence+Alpha) Error opening " << fnAuxBCR << "." << std::endl;
        exit (EXIT_FAILURE);
    }
   
    assert (fread(&lengthTot_plus_eof,sizeof(dataTypeNChar),1,InfoFileBCR) ==1);
    assert (fread(&nText,sizeof(dataTypeNSeq),1,InfoFileBCR) ==1);
    assert (fread(&sizeAlpha,sizeof(dataTypedimAlpha),1,InfoFileBCR) ==1);
    assert (fread(&lengthRead,sizeof(dataTypelenSeq),1,InfoFileBCR) ==1);

    
    uchar *alphaBCR = new uchar[sizeAlpha];
    alphaInverse = new dataTypedimAlpha[sizeAlpha];

    //set alpha and alphaInverse
    cout << "Alpha and code:\n";
    for (dataTypedimAlpha i = 0; i < sizeAlpha; i++) {
        assert ( fread(&alphaBCR[i],sizeof(uchar),1,InfoFileBCR) == 1);
        alpha[alphaBCR[i]] = i;
        alphaInverse[i]=alphaBCR[i];
        cout << (int)alpha[alphaBCR[i]] <<  "\t" <<  (int)alphaBCR[i] << "\t" <<  (int)alphaInverse[i]<< "\n";
    }
    fclose(InfoFileBCR);
    delete [] alphaBCR;
    
    //set tableOcc
    tableOcc = new dataTypeNChar*[sizeAlpha];
    //Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        tableOcc[j] = new dataTypeNChar[sizeAlpha];
    }
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++)
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
            tableOcc[j][h]=0;
        
    string fnAuxTable = filenameInfo + ".table\0";
    FILE* InFileTableBCR = fopen(fnAuxTable.c_str(), "rb");
    if (InFileTableBCR==NULL) {
        std::cerr << "(reading table occ) Error opening " << fnAuxTable << std::endl;
        exit (EXIT_FAILURE);
    }
    
    //set freq
    for (dataTypedimAlpha j = 0 ; j < SIZE_ALPHA-1; j++) {
        freq[j]=0;
    }
    freq[SIZE_ALPHA-1]=0;

    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        //Update tableOcc with the symbols of the previous partial BWT files
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++) {
            assert (fread(&tableOcc[j][h],sizeof(dataTypeNChar),1,InFileTableBCR) == 1);
            freq[(unsigned int)alphaInverse[j]] += tableOcc[j][h];
        }
    }
    fclose(InFileTableBCR);

    std::cout << "\nFrom " << fnAuxBCR << " file:\n";
    std::cout << "\tNumber of sequences: " << nText << "\n";
    std::cout << "\tTotal length (with $): " << lengthTot_plus_eof << "\n";
    std::cout << "\tSize alpha: " << (int)sizeAlpha << "\n";
    std::cout << "\tMaximum Length of the sequences: " << (int)lengthRead << "\n";
        
    //#if DEBUG == 1
    std::cout << "\nFrom " << fnAuxTable << " file (TableOcc and freq):\n";
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
            std::cout << tableOcc[j][h] << "\t";
        std::cout << "freq[" << alphaInverse[j] << "]= " << freq[(unsigned int)alphaInverse[j]] <<  " \n";
    }
        //#endif

    return 1;
}
    

    
int BCRdecode::computeVectorUnbuildBCR()
    {
        numBlocksInPartialBWT.resize(sizeAlpha);
        //Set number of blocks for each BWT-partial
        for (dataTypedimAlpha x = 0 ; x <= sizeAlpha-1; x++) {
            numBlocksInPartialBWT[x] = (dataTypeNChar)ceil((long double)freq[alphaInverse[x]]/DIMBLOCK);
            //#if DEBUG==1
                std::cerr << "computeVectorUnbuildBCR: " << "freq[" << (unsigned int)x << "]= " << freq[alphaInverse[x]] << " and numBlocksInPartialBWT[" << (unsigned int)x << "]= " << numBlocksInPartialBWT[x] << "\n";
            //#endif
        }


        // Start by allocating an array for array of arrays
        
        vectorOcc.resize(sizeAlpha);    //For each BWT-partial

        #if OMP==0
            time_t start,end;
            double dif;
        #endif

        
        // alphaInverse[x] is the symbol to which correspond bwt_x

        //For each BWT-partial
        //Read BWT-partials in parallel
        
        #if OMP
            double d_total = omp_get_wtime();
        
            vector < FILE * > InFileBWT;
            InFileBWT.resize(numthreads);
        
        #else
            static FILE *InFileBWT;
        #endif
        
        dataTypedimAlpha x=0;
        
        #if OMP
            #pragma omp parallel for default(shared) private(x) num_threads(numthreads) schedule(dynamic, 1)
        #endif
            for ( x=0; x < sizeAlpha; x++)
            {
                #if OMP
                    int t = omp_get_thread_num(); //t is the thread ID
                
                    InFileBWT[t] = openFilePartialIn(x);
                    fseek(InFileBWT[t], 0, SEEK_SET);
                #else
                    //x is the cycle index
                    InFileBWT = openFilePartialIn(x);
                    fseek(InFileBWT, 0, SEEK_SET);
                #endif
                
                // Allocate an array for each block of BWT-partial
                vectorOcc[x].resize(sizeAlpha);
                
                // Allocate an array of integers for each element of vectorOcc[x]
                for (dataTypedimAlpha y = 0 ; y < sizeAlpha; y++)   {      //For each block
                    vectorOcc[x][y].resize(numBlocksInPartialBWT[x],0);
                }
                
                uchar *bufBlock = new uchar[DIMBLOCK];

                dataTypeNChar numBlock = 0;
                dataTypeNChar num_read = 0;

                
                #if OMP
                    double tr_start=omp_get_wtime();
                #else
                    time (&start);
                #endif
                
                #if OMP==0 && DEBUG==1
                    std::cerr << "\ncomputeVectorUnbuildBCR: " << " start= " << start << "= " << " and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
                #endif

                //Read DIMBLOCK symbols in BWT-partial
                #if OMP
                    while( ( (num_read =  readOnFilePartial(bufBlock, DIMBLOCK, InFileBWT[t]) ) && (num_read > 0) )  &&  (numBlock < numBlocksInPartialBWT[x]))   //Added check on numBlocks
                #else
                    while( ( (num_read =  readOnFilePartial(bufBlock, DIMBLOCK, InFileBWT) ) && (num_read > 0) )  &&  (numBlock < numBlocksInPartialBWT[x]))   //Added check on numBlocks
                #endif
                {
                    for (dataTypeNChar i=0; i<num_read; i++) {
                        vectorOcc[x][alpha[(unsigned int)(bufBlock[i])]][numBlock]++;
                        //std::cerr << "---x = " << (unsigned int)x << " alpha " << (unsigned int)(bufBlock[i]) << " -numBlock " << numBlock << " vectorOcc is " << vectorOcc[x][alpha[(unsigned int)(bufBlock[i])]][numBlock] << ".\n";
                    }
                    numBlock++;
                    
                }//end-while
                
                #if OMP
                if ( !feof(InFileBWT[t]) && (numBlock > numBlocksInPartialBWT[x])) {
                #else
                if ( !feof(InFileBWT) && (numBlock > numBlocksInPartialBWT[x])) {
                #endif
                    std::cerr << "computeVectorUnbuildBCR: Error - The file contains more blocks than allocates.\n" ;
                    exit(1);
                }
                
                //For FIXED x
                //Compute the sum cumulative for each BWT-partial
                for (dataTypedimAlpha z = 0 ; z < sizeAlpha; z++)  {      //For each symbol z
                    for(dataTypeNChar y = 1; y < numBlocksInPartialBWT[x] ; y++)   {      //For each block y>1 of partial-BWT x
                        vectorOcc[x][z][y]=vectorOcc[x][z][y-1] + vectorOcc[x][z][y];   //Sum the previous one: ie Blcok y and block y-1
                    }
                }
            
                    //VectorOcc[x] stores info similar to row tableOcc[x] but keeping track of symbol occurrences into blocks
                    //N.B. VectorOcc[x][z][y] stores the total number of z-occurrences in the BWT-partial corresponding to alpha[x] symbol up to the y-th block
            
            
                #if OMP
                #pragma omp critical
                    {
                        std::cerr << "TIME THREAD " << t << " = " << omp_get_wtime()- tr_start << "(in seconds) and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
                    }
                #else
                    time (&end);
                    dif = difftime (end,start);
                    
                    std::cerr << "computeVectorUnbuildBCR: the cycle for " << "symbol = " << (int)x << " tooks " << dif << " seconds" << " and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
                #endif
                    
                #if OMP==0 && DEBUG==1
                    std::cerr << "computeVectorUnbuildBCR: " << " end= " << end << "]= " << " and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
                #endif
                
                delete [] bufBlock;
                    
                #if OMP
                closeFilePartial(InFileBWT[t]);
                #else
                closeFilePartial(InFileBWT);
                #endif
                
        }//end-for

        #if OMP
                std::cout << "Time: " << omp_get_wtime()-d_total << std::endl;
        #endif
        
        #if DEBUG==1
            for (dataTypedimAlpha x = 0 ; x < sizeAlpha; x++) {
                std::cerr << "x = " << (unsigned int)x << " For the " << alphaInverse[x] << "-BWT-partial: the #symbols is " << freq[alphaInverse[x]] << ".\n";
                std::cerr << "Number of blocks of the symbol " << numBlocksInPartialBWT[x] << "\n";
                for(dataTypedimAlpha z = 0; z < sizeAlpha; ++z) {
                    std::cerr << "Symbol " << (unsigned int)z << ":\t";
                    for(dataTypeNChar y = 0; y < numBlocksInPartialBWT[x]; ++y) {
                        std::cerr << (int)vectorOcc[x][z][y];
                    }
                    std::cerr << "\n";
                }
                std::cerr << "\n";
            }
        #endif
        
        return 1;
    }

int BCRdecode::findBlockToRead(dataTypeNChar *counters, dataTypedimAlpha currentPile, dataTypeNChar *toRead, dataTypeNChar *numBlock) {
        //Find the block numblock, where the position toRead is
    
        *numBlock = (dataTypeNChar)floor((long double)((*toRead-1)/DIMBLOCK)) ;  //The smallest integral value NOT less than x.
        //if (*numBlock >= numBlocksInPartialBWT[currentPile])
        //std::cerr << "     findBlockToRead: numBlock " << *numBlock << " and numBlocksInPartialBWT["<<(unsigned int)currentPile<<"]= " << numBlocksInPartialBWT[currentPile] << "\n";
        assert(*numBlock < numBlocksInPartialBWT[currentPile]);
    
        if (*numBlock > 0) {
           for (dataTypedimAlpha r=0; r<sizeAlpha; r++)
              counters[r] =  vectorOcc[currentPile][r][(*numBlock)-1];   //vectorOcc is indexed by 0, so we have numBlock-1
           *toRead = *toRead - (*numBlock*DIMBLOCK);  //Number of symbols that we must read yet. it could be = DIMBLOCK
        }
        return 1;
    }
    
    
int BCRdecode::buildFreq(string fileWithExt) {
   
    //Open LCP and DA and BWT files
    string fnBWT = string(fileWithExt) + ".ebwt";
    std::cerr << "Split " << fnBWT << " file in BWT/LCP/DA/QS partial file." << std::endl;
    FILE *InBWT = fopen(fnBWT.c_str(), "rb");
    if (InBWT==NULL) {
        std::cerr << "Error opening " << fnBWT << "!" << std::endl;
        exit (EXIT_FAILURE);
    }
    fseek(InBWT, 0, SEEK_SET);
    #if BUILD_LCP == 1
        string fnLCP = string(fileWithExt) +".lcp\0";
        FILE *InLCP = fopen(fnLCP.c_str(), "rb");
        if (InLCP==NULL) {
            std::cerr << "Error opening " << fnLCP << "." << std::endl;
            exit (EXIT_FAILURE);
        }
        fseek(InLCP, 0, SEEK_SET);
    #endif
    #if BUILD_DA == 1
        string fnDA = string(fileWithExt)+".da\0";
        FILE *InDA = fopen(fnDA.c_str(), "rb");
        if (InDA==NULL) {
            std::cerr << "Error opening " << fnDA << "." << std::endl;
            exit (EXIT_FAILURE);
        }
        fseek(InDA,0, SEEK_SET);
    #endif
    #if USE_QS == 1
        string fnQS = string(fnBWT) + ".qs";
        FILE *InQS = fopen(fnQS.c_str(), "rb");
        if (InQS==NULL) {
            std::cerr << "Error opening " << fnQS << "." << std::endl;
            exit (EXIT_FAILURE);
        }
        fseek(InQS,0, SEEK_SET);
    #endif
    
    dataTypeNChar numEle=0;
    
    for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
        freq[z]=0;
    freq[SIZE_ALPHA-1]=0;
    
    
    //First reading in order to find the alphabet
    std::cerr << "Find the alphabet by reading " << fnBWT << " file" << std::endl;
    uchar bwt;
    //uchar c;
    dataTypeNChar numcharBWT;
    numcharBWT = fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT);
    
    #if BUILD_LCP == 1
        dataTypelenSeq lcp;
        dataTypeNChar numcharLCP;
        numcharLCP = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
        assert(numcharLCP==numcharBWT);
        dataTypelenSeq minLCP=lcp, maxLCP=lcp;
    #endif
    #if BUILD_DA == 1
        dataTypeNChar numcharDA;
        dataTypeNSeq idSeq;
        numcharDA = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
        assert(numcharBWT==numcharDA);
    #endif
    #if USE_QS == 1
        dataTypeNChar numcharQS;
        uchar qs;
        numcharQS = fread(&qs,sizeof(uchar),1,InQS);
        assert(numcharBWT==numcharQS);
    #endif
    
    
    freq[(unsigned int)(bwt)]=1;
    numEle++;
    
    //set freq
    dataTypeNSeq numSeq=0;
    numcharBWT=1;
    while ( ( numcharBWT = fread(&bwt,sizeof(dataTypedimAlpha),1, InBWT) ) && ( numcharBWT > 0 )  ) {
        #if BUILD_LCP == 1
            numcharLCP = fread(&lcp,sizeof(dataTypelenSeq),1,InLCP);
            if (minLCP > lcp)
                    minLCP=lcp;
            if (maxLCP < lcp)
                    maxLCP = lcp;
            assert(numcharBWT == numcharLCP);
        #endif
        #if BUILD_DA == 1
            numcharDA = fread(&idSeq,sizeof(dataTypeNSeq),1,InDA);
            assert(numcharBWT == numcharDA);
        #endif
        #if USE_QS == 1
            numcharQS = fread(&qs,sizeof(uchar),1,InQS);
            assert(numcharBWT == numcharQS);
        #endif
        
        if (bwt == TERMINATE_CHAR)
            numSeq++;

        freq[(unsigned int)(bwt)]++;
        
        numEle++;
        //cerr << numEle << " " << bwt << " " << (int)qs << "\n";
    }//end-while
    
    //set nText
    nText = numSeq;
    
    //set alpha and alphaInverse
    alphaInverse = new dataTypedimAlpha[SIZE_ALPHA];
    sizeAlpha=0;
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
        if (freq[i] > 0) {
            alpha[i] = sizeAlpha;
            alphaInverse[sizeAlpha]=i;
            sizeAlpha++;
        }
    if (freq[SIZE_ALPHA-1] > 0) {
        alpha[SIZE_ALPHA-1] = sizeAlpha;
        alphaInverse[sizeAlpha]=SIZE_ALPHA-1;
        sizeAlpha++;
    }
    
    std::cerr << "\nFrom .ebwt file:\n";
    std::cerr << "\tNumber of sequences: " << numSeq << "\n";
    std::cerr << "\tNumber of symbols in the input file: " << numEle << "\n";
    std::cerr << "\tSize alpha: " << (int)sizeAlpha << "\n";
    #if BUILD_LCP == 1
        std::cerr << "\tminLCP: " << minLCP << "\n";
        std::cerr << "\tmaxLCP: " << maxLCP << "\n";
    #endif

    fclose(InBWT);
    #if BUILD_LCP == 1
        fclose(InLCP);
    #endif
    #if BUILD_DA == 1
        fclose(InDA);
    #endif
    #if USE_QS == 1
        fclose(InQS);
    #endif
    
    return 1;
}
    
int BCRdecode::splitIntoPartial(string fileWithExt, int mode) {
    
    //Open BWT, LCP, DA, QS files
    string fnBWT = string(fileWithExt) + ".ebwt";
    vector < FILE * > InBWT;
    
    int t=0;
    
    #if OMP
    InBWT.resize(numthreads);
    for(t=0;t<numthreads; t++)
    #else
    InBWT.resize(1);
    #endif
    {
        InBWT[t] = fopen(fnBWT.c_str(), "rb");
        if (InBWT[t]==NULL) {
            std::cerr << "Error opening " << fnBWT << "!" << std::endl;
            exit (EXIT_FAILURE);
        }
    }

    #if BUILD_LCP == 1
        string fnLCP = string(fileWithExt) +".lcp\0";
        vector < FILE * > InLCP;
        #if OMP
        InLCP.resize(numthreads);
        for(t=0;t<numthreads; t++)
        #else
        InLCP.resize(1);
        #endif
        {
            InLCP[t] = fopen(fnLCP.c_str(), "rb");
            if (InLCP[t]==NULL) {
                std::cerr << "Error opening " << fnLCP << "." << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    #endif
    
    #if BUILD_DA == 1
        string fnDA = string(fileWithExt)+".da\0";
        vector < FILE * > InDA;
        #if OMP
        InDA.resize(numthreads);
        for(t=0;t<numthreads; t++)
        #else
        InDA.resize(1);
        #endif
        {
            InDA[t] = fopen(fnDA.c_str(), "rb");
            if (InDA[t]==NULL) {
                std::cerr << "Error opening " << fnDA << "." << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    #endif
    
    #if USE_QS == 1
        string fnQS = string(fnBWT) + ".qs";
        vector < FILE * > InQS;
        #if OMP
        InQS.resize(numthreads);
        for(t=0;t<numthreads; t++)
        #else
        InQS.resize(1);
        #endif
        {
            InQS[t] = fopen(fnQS.c_str(), "rb");
            if (InQS[t]==NULL) {
                std::cerr << "Error opening " << fnQS << "." << std::endl;
                exit (EXIT_FAILURE);}
        }
    #endif
    
    if (mode == 1) {    //allocate tableOcc
        tableOcc = new dataTypeNChar*[sizeAlpha];
        //Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
        for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
            tableOcc[j] = new dataTypeNChar[sizeAlpha];
        }
        for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++)
            for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
                tableOcc[j][h]=0;
    }
    
    //BCR partial files
    std::cerr << "Build the partial BCR files for eBWT/LCP/DA/QS." << std::endl;

    vector < FILE *> OutFileBWT;
    #if OMP
        OutFileBWT.resize(numthreads);
    #else
        OutFileBWT.resize(1);
    #endif
        

    #if BUILD_LCP == 1
        vector < FILE *> OutFileLCP;
        #if OMP
        OutFileLCP.resize(numthreads);
        #else
        OutFileLCP.resize(1);
        #endif
    #endif
        
    #if BUILD_DA == 1
        vector < FILE *> OutFileDA;
        #if OMP
        OutFileDA.resize(numthreads);
        #else
        OutFileDA.resize(1);
        #endif
    #endif
        
    #if USE_QS == 1
        vector < FILE *> OutFileQS;
        #if OMP
        OutFileQS.resize(numthreads);
        #else
        OutFileQS.resize(1);
        #endif
    #endif
        
    dataTypeNChar numEle=0;
    
    //sumFreq si to split ebwt
    vector<dataTypeNChar> sumFreq;
    sumFreq.resize(sizeAlpha+1);
    
    sumFreq[0]=0;
    for (dataTypedimAlpha j = 1 ; j <= sizeAlpha; j++)
        sumFreq[j]=sumFreq[j-1]+freq[alphaInverse[j-1]];
    
    size_t currentPile;
    #if OMP
#pragma omp parallel for default(shared) private(currentPile) firstprivate(t) num_threads(numthreads) schedule(dynamic, 1) reduction(+:numEle)
    #endif
        for (currentPile = 0 ; currentPile < sizeAlpha; ++currentPile) {
            
            assert(freq[alphaInverse[currentPile]] > 0);

            #if OMP
                t = omp_get_thread_num();//id_thread
                double start = omp_get_wtime();
            #endif
            
            //Open BCR partial files
            dataTypeNChar numcharBWT, numWrite;
            uchar *bufferBWT = new uchar[DIMBLOCK];
            
            string fnOutBWT = "bwt_" + to_string(currentPile)+ext;
            OutFileBWT[t] = fopen(fnOutBWT.c_str(), "wb");
            if (OutFileBWT[t]==NULL)
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "Error opening: " << fnOutBWT << std::endl;
                exit (EXIT_FAILURE);
            }
            fseek(InBWT[t], sumFreq[currentPile]*sizeof(uchar), SEEK_SET);
            
            #if BUILD_LCP == 1
            string fnOutLCP = "lcp_" + to_string(currentPile)+ext;
            OutFileLCP[t] = fopen(fnOutLCP.c_str(), "wb");
            if (OutFileLCP[t]==NULL)
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "Error opening: " << fnOutLCP << std::endl;
                exit (EXIT_FAILURE);
            }
            fseek(InLCP[t], sumFreq[currentPile]*sizeof(dataTypelenSeq), SEEK_SET);
            dataTypeNChar numcharLCP;
            dataTypelenSeq *bufferLCP = new dataTypelenSeq[DIMBLOCK];
            #endif
            
            #if BUILD_DA == 1
            string fnOutDA = "da_" + to_string(currentPile)+ext;
            OutFileDA[t] = fopen(fnOutDA.c_str(), "wb");
            if (OutFileDA[t] ==NULL)
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "Error opening: " << fnOutDA << std::endl;
                exit (EXIT_FAILURE);
            }
            fseek(InDA[t], sumFreq[currentPile]*sizeof(dataTypeNSeq), SEEK_SET);
            dataTypeNChar numcharDA;
            dataTypeNSeq *bufferDA = new dataTypeNSeq[DIMBLOCK];
            #endif
            
            #if USE_QS == 1
            string fnOutQS = "bwt_qs_" + to_string(currentPile)+ext;
            OutFileQS[t] = fopen(fnOutQS.c_str(), "wb");
            if (OutFileQS[t] ==NULL)
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "Error opening: " << fnOutQS << std::endl;
                exit (EXIT_FAILURE);
            }
            fseek(InQS[t], sumFreq[currentPile]*sizeof(uchar), SEEK_SET);
            dataTypeNChar numcharQS;
            uchar *bufferBWTQS = new uchar[DIMBLOCK];
            #endif

            //Build BCR partial files
            dataTypeNChar j = sumFreq[currentPile];
            
            dataTypeNChar toRead = DIMBLOCK;  //read DIMBLOCK symbols per time
            
            while (j < sumFreq[currentPile+1]) {
                
                if( (sumFreq[currentPile+1] - j) < toRead) //check remaining symbols to read
                    toRead = sumFreq[currentPile+1] - j;
                
                //read toRead symbols
                numcharBWT = fread(bufferBWT, sizeof(uchar),toRead, InBWT[t]);
                //write
                numWrite = fwrite(bufferBWT, sizeof(uchar), toRead, OutFileBWT[t]);
                assert(numcharBWT == numWrite);
                
                //set tableOcc
                if(mode == 1){
                    //counting the number of occurrences in BWT of the currentPile
                    uchar bwt;
                    for(dataTypeNChar i = 0; i < numcharBWT; i++ ){
                        bwt = bufferBWT[i];
                        tableOcc[(unsigned int)currentPile][alpha[(unsigned int)bwt]]++;
                    }
                }

                
                #if BUILD_LCP == 1
                numcharLCP = fread(bufferLCP, sizeof(dataTypelenSeq),toRead, InLCP[t]);
                assert(numcharLCP==numcharBWT);
                numWrite = fwrite (bufferLCP, sizeof(dataTypelenSeq), toRead , OutFileLCP[t]);
                assert(numcharLCP == numWrite);
                #endif
                
                #if BUILD_DA == 1
                numcharDA = fread(bufferDA, sizeof(dataTypeNSeq),toRead, InDA[t]);
                
                if (numcharDA!=numcharBWT)
                    std::cerr <<  "numcharDA "<< numcharDA << " and numcharBWT" << numcharBWT << "\n";
                
                assert(numcharDA==numcharBWT);
                numWrite = fwrite (bufferDA, sizeof(dataTypeNSeq),toRead, OutFileDA[t]);
                assert(numcharDA == numWrite);
                #endif
                
                #if USE_QS == 1
                numcharQS = fread(bufferBWTQS, sizeof(uchar),toRead, InQS[t]);
                
                if (numcharQS!=numcharBWT)
                    std::cerr <<  "numcharQS "<< numcharQS << " and numcharBWT" << numcharBWT << "\n";
                
                assert(numcharQS==numcharBWT);
                numWrite = fwrite (bufferBWTQS, sizeof(uchar),toRead, OutFileQS[t]);
                assert(numcharQS == numWrite);
                #endif
                
                numEle++;
                
                j+=toRead;
                
            }  //end-for
            
            //Close partial files
            fclose(OutFileBWT[t]);
            #if BUILD_LCP == 1
                fclose(OutFileLCP[t]);
            #endif
            #if BUILD_DA == 1
                fclose(OutFileDA[t]);
            #endif
            #if USE_QS == 1
                fclose(OutFileQS[t]);
            #endif
            
            #if OMP
            #pragma omp critical
            {
                std::cerr << "splitIntoPartial: THREAD = " << t << " tooks " << omp_get_wtime()-start << " seconds " << "\n\n";
            }
            #endif
            
        }  //end-for
    
    //Close BWT, LCP, DA, QS files
    #if OMP
    for(t=0;t<numthreads; t++)
    #endif
    {
        fclose(InBWT[t]);
        #if BUILD_LCP == 1
            fclose(InLCP[t]);
        #endif
        #if BUILD_DA == 1
            fclose(InDA[t]);
        #endif
        #if USE_QS == 1
            fclose(InQS[t]);
        #endif
    }
        
        
    std::cerr <<  "The total number of elements is " << numEle << "\n";
        
    return 1;
    
}
        
        
    
    //Multiple Decoding the sequences (Build reverse sequence)
    //Reconstruct m sequences backwards by threading through the mapping and reading the characters off of L.
    //fileInput is the input file
    //fileOutBWT is the suffix of the filename of the partial BWTs
    //fileOutCyc is the prefix of the lengthRead-filename (transpose texts: cyc.i.txt)
        
int BCRdecode::decodeBCRmultipleReverse(string fileInput)
    {
        vectInsTexts.resize(nText,0);
        
        sortElement pair;
        
        #if OMP
            int j;
            vectPairParallel.resize(numthreads);
            vectPairParallelTmp.resize(numthreads);
            for (j=0; j<numthreads; j++){
                vectPairParallel[j].resize(sizeAlpha);
                vectPairParallelTmp[j].resize(sizeAlpha);
            }
            //vectPairParallel[j][pile] contains the pairs for the j-th thread
        #else
            vectPair.resize(sizeAlpha);
            vectPairTmp.resize(sizeAlpha);
        #endif

        //Position of $ in F
        dataTypedimAlpha pile=alpha[(unsigned int)(TERMINATE_CHAR)];
        
        #if OMP==0
            for (dataTypeNSeq g = 0 ; g < nText; g++) {
                pair.posN = g + 1;
                pair.seqN = g;
                vectPair[pile].push_back(pair);
            }
            #if DEBUG==1
                std::cerr << "The Initial triples of " << (int)pile << " in first column are!"<< std::endl;
                for (std::list<sortElement>::const_iterator it = vectPair[pile].begin(); it != vectPair[pile].end(); ++it)
                    std::cout <<  "\t" <<  it->posN << "\t" << it->seqN << std::endl;
                std::cerr << "vectSizeCurrentPile[" << (int)pile << "]= " << vectSizeCurrentPile[pile] << std::endl;
            #endif
        #endif

        //newSymb
        uchar *newSymb = NULL;

        #if USE_QS == 1
            newSymb = new uchar[2*nText];
        #else
            newSymb = new uchar[nText];
        #endif
        
        #if OMP==0
            static FILE *InfileOutDecodeCyc;
            string filename;
            time_t start,end;
        #endif
        
        #if OMP
            vector < FILE * > InfileOutDecodeCyc;
            InfileOutDecodeCyc.resize(numthreads);
        
            //numPairsIn is # pairs in vectPairParallel
            dataTypeNChar numPairsIn=(dataTypeNChar)ceil((long double)nText/numthreads);
            cout << "numPairsIn:" << numPairsIn << endl;
        
            #pragma omp parallel default(shared) private(pair) firstprivate(numPairsIn) num_threads(numthreads)
            {
                int tid=omp_get_thread_num();//id_thread
                double start=omp_get_wtime();
                
                dataTypeNSeq s_ind = tid*numPairsIn; //start index
                
                #pragma omp master
                {
                    std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
                }
                if(tid==numthreads-1){//last thread
                    numPairsIn=nText-s_ind; //set numPairsIn for last thread
                }
                //push pairs vectPairParallel[tid]
                for (dataTypeNSeq ind=s_ind; ind < s_ind+numPairsIn; ind++) {
                    pair.posN = ind + 1;
                    pair.seqN = ind;
                    vectPairParallel[tid][pile].push_back(pair);
                }
                string filename;
                #if USE_QS == 1
                    string filenameQS;
                #endif
                
                dataTypelenSeq m;
                //As we recover the symbol in reverse order, I store the first found symbol in cyc.(length-1).txt file
                //and the last found symbol in cyc.0.txt file
                for (m = lengthRead ; m > 0 ; m--) {

                    assert ( RecoverNsymbolsReverseByVector(tid, fileInput, newSymb) == 1);

                    filename = fileOutCyc + to_string(tid) + "_" + to_string(m-1) +".txt";
                    //open file
                    InfileOutDecodeCyc[tid] = fopen(filename.c_str(), "wb");
                    if (InfileOutDecodeCyc[tid]==NULL) {
                        #pragma omp critical
                        {
                            std::cerr << "decodeBCRmultipleReverse: could not open file " << filename << " !" << std::endl;
                            exit (EXIT_FAILURE);
                        }
                    }

                    #if USE_QS == 1
                        dataTypeNChar numcharWrite = fwrite (&newSymb[2*s_ind], sizeof(uchar), 2*numPairsIn , InfileOutDecodeCyc[tid]);
                        assert( numcharWrite == 2*numPairsIn); // we should always read the same number of characters
                    #else
                        dataTypeNChar numcharWrite = fwrite (&newSymb[s_ind], sizeof(uchar), numPairsIn , InfileOutDecodeCyc[tid]);
                        assert( numcharWrite == numPairsIn); // we should always read the same number of characters
                    #endif
                    
                    //close file
                    fclose(InfileOutDecodeCyc[tid]);
                    
                }//end-for
                
                #pragma omp critical
                {
                    std::cerr << "decodeBCRmultipleReverse: THREAD = " << tid << " tooks " << omp_get_wtime()-start << " seconds " << "\n\n";
                }
                
            }//end-pragma
            #else
                for (dataTypelenSeq m = lengthRead ; m > 0 ; m--) {
                    
                    time (&start);
                    
                    assert ( RecoverNsymbolsReverseByVector(fileInput, newSymb) == 1);
                    
                    filename = fileOutCyc + to_string(m-1) +".txt";
                    //sprintf (filename, "%s%u.txt", fileOutCyc, m-1);
                    InfileOutDecodeCyc = fopen(filename.c_str(), "wb");
                    if (InfileOutDecodeCyc==NULL) {
                        std::cerr << "decodeBCRmultipleReverse: could not open file " << filename << " !" << std::endl;
                        exit (EXIT_FAILURE);
                    }

                    #if USE_QS == 1
                        dataTypeNChar numcharWrite = fwrite (newSymb, sizeof(uchar), 2*nText , InfileOutDecodeCyc);
                        assert( numcharWrite == 2*nText); // we should always read the same number of characters
                    #else
                        dataTypeNChar numcharWrite = fwrite (newSymb, sizeof(uchar), nText , InfileOutDecodeCyc);
                        assert( numcharWrite == nText); // we should always read the same number of characters
                    #endif
                    
                    fclose(InfileOutDecodeCyc);
                    
                    time (&end);
                    #if DEBUG
                        std::cerr << "decodeBCRmultipleReverse: the cycle for position = " << (int)m << " tooks " << difftime (end,start) << " seconds" << "\n\n";
                    #endif
                    
                }//end-for
            #endif
        
        delete [] newSymb;

        return 1;
    }
    
    
//It is used to reconstruct m sequences backwards by threading through the mapping and reading the characters off of L.
#if OMP
        int BCRdecode::RecoverNsymbolsReverseByVector(int t_id, string fileInput,  uchar * newSymb)    {
#else
        int BCRdecode::RecoverNsymbolsReverseByVector(string fileInput,  uchar * newSymb)    {
#endif
    
    string filename;
    FILE *InFileBWT;
    #if USE_QS
        FILE *InFileQS=NULL;
        string filenameQS;
    #endif
    
    dataTypeNChar toRead = 0;
    dataTypeNChar *counters = new dataTypeNChar[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    uchar *bufferBlock = new uchar[DIMBLOCK];
            
    dataTypedimAlpha currentPile=0, newCurrentPile=0;
    
    dataTypeNChar posInPile=0;
    
    sortElement pair, newPair;

    while(currentPile < sizeAlpha){
        
        //Set counters to 0
        for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
            counters[i]=0;
            
        filename = fileOutBwt + to_string((int)(currentPile)) + ext;
        
        //Open BWT-partial
        InFileBWT = fopen(filename.c_str(), "rb");
        if (InFileBWT==NULL) {
            #if OMP
                #pragma omp critical
            #endif
            {
                std::cerr << "RecoverNsymbolsReverseByVector: " << filename << " file error opening " << std::endl;
                exit (EXIT_FAILURE);
            }
        }
            
        #if USE_QS
            filenameQS = fileOutBwt + "qs_" + to_string((int)(currentPile)) + ext;
            InFileQS = fopen(filenameQS.c_str(), "rb");
            if (InFileQS==NULL) {
                #if OMP
                    #pragma omp critical
                #endif
                {
                    std::cerr << "RecoverNsymbolsReverseByVector: " << filenameQS << " file error opening " << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
        #endif

        uchar foundSymbol, foundQual='\0';
        
        dataTypeNChar numBlock = 0;
           
        posInPile=0; //posInPile is the position in vectPair[currentPile]
        
        #if OMP
        
            #if DEBUG
                #pragma omp critical
                {
                    cout << "t_id: " << t_id << ", currentPile: " << (int)currentPile << ", vectPairParallel[t_id][currentPile].size: " << vectPairParallel[t_id][currentPile].size() << endl;
                }
            #endif
        
        while (!vectPairParallel[t_id][currentPile].empty() ) {
            //Pick up a pair
            pair = vectPairParallel[t_id][currentPile].front();
            //Remove that pair
            vectPairParallel[t_id][currentPile].pop_front();
        #else
        while (!vectPair[currentPile].empty() ) {
            //Pick up a pair
            pair = vectPair[currentPile].front();
            //Remove that pair
            vectPair[currentPile].pop_front();
        #endif

            if (vectInsTexts[pair.seqN]==0) { //Retrieve the BWT-symbol to insert in newSymb[i]
                    
                #if DEBUG==1 && OMP==0
                    std::cerr << "Inside While Sequence number k= " << posInPile << "\n";
                    std::cerr << "\t" << " P["<< posInPile <<"]=" << (dataTypeNChar)pair.posN << " N["<< posInPile <<"]=" << (dataTypeNSeq)pair.seqN << "\n";
                #endif
                    
                //For any character (of differents sequences) in the same pile
                foundSymbol = '\0';
                #if USE_QS
                    foundQual = '\0';
                #endif
                
                toRead = pair.posN; //skip pair.posN symbols

                //set to zero counters
                for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                        counters[i]=0;
                    
                if (toRead > 0) {
                        //if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks preceding the block where the position toRead is.
                    assert (findBlockToRead(counters, currentPile, &toRead, &numBlock) == 1);
                }

                if (toRead > 0 && toRead <= DIMBLOCK) {   //If toRead == DIMBLOCK, because I can need to known foundSymbol character
                       
                    fseek (InFileBWT, numBlock*DIMBLOCK, 0);
                    #if USE_QS
                        fseek (InFileQS, numBlock * DIMBLOCK, 0 );
                    #endif
                        
                    #if DEBUG==1 && OMP==0
                        std::cerr << "RecoverNsymbolsReverseByVector Move file to the position " << numBlock*DIMBLOCK <<  "\n";
                    #endif
                    
                    //Read BWT-symbol
                    dataTypeNChar numchar;
                    
                    #if ( USE_QS == 1 )
                        dataTypeNChar toReadQS = toRead;
                    #endif
                    
                    numchar = fread(bufferBlock,sizeof(uchar),toRead,InFileBWT);
                    assert(numchar == toRead);  // we should always read/write the same number of characters
                    
                    #if DEBUG == 1 && OMP==0
                        std::cerr << "toRead " << toRead << " is less than DIMBLOCK " << DIMBLOCK << "\n";
                        std::cerr << "numchar (should be equal to toRead): " << numchar << "\n";
                    #endif
                    
                    #if DEBUG && OMP
                    #pragma omp critical
                    {
                        cout << "t_id: " << t_id << ", numchar: " << numchar << ", toRead: " << toRead << endl;
                    }
                    #endif
                    
                    foundSymbol = bufferBlock[numchar-1];//The symbol is in the last position in the partial BWT that we have read.
                    
                    #if DEBUG == 1 && OMP==0
                        std::cerr << "** counters before:\t";
                        for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                            std::cerr << " " << counters[i];
                        std::cerr << "\n";
                    #endif
                    
                    //For each symbol in the buffer, it updates the number of occurrences into counters
                    for (dataTypeNChar r=0; r<numchar; r++)
                        counters[alpha[(unsigned int)bufferBlock[r]]]++;//increment the number of letter symbol into counters
                    
                    #if DEBUG == 1 && OMP==0
                        std::cerr << "counters after:\t";
                        for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                            std::cerr << " " << counters[i];
                        std::cerr << "\n**\n";
                    #endif
                    
                    #if ( USE_QS == 1 )
                        if (toReadQS > 1) {
                            fseek(InFileQS, (toReadQS - 1), SEEK_CUR );
                        }
                        assert (fread( &foundQual, sizeof( uchar ), 1, InFileQS ) == 1);
                        
                        size_t pos1 = ftell( InFileBWT );
                        size_t pos2 = ftell( InFileQS );
                        assert( pos1 == pos2 );
                    #endif

                    #if DEBUG==1 && OMP==0
                        std::cerr << "RecoverNsymbolsReverseByVector foundSymbol " << (unsigned int)foundSymbol << " " << (char)foundSymbol <<  "\n";
                        #if ( USE_QS == 1 )
                            std::cerr << "RecoverNsymbolsReverseByVector foundQual " << (unsigned int)foundQual <<  "\n";
                        #endif
                    #endif
                        
                }  //end-if (toRead <= DIMBLOCK)
                    
                #if DEBUG==1 && OMP==0
                    std::cerr << "counters after FirstVector:\t";
                    for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                            std::cerr << " " << counters[i];
                    std::cerr << "\n";
                    std::cerr << "pair.seqN = " << pair.seqN << " Symbol = " << foundSymbol << "\n";
                #endif
  
                //write foundSymbol and foundQual
                #if USE_QS
                    newSymb[2*pair.seqN] = foundSymbol;
                    newSymb[2*pair.seqN+1] = foundQual;
                #else
                    newSymb[pair.seqN] = foundSymbol;
                #endif
                
                //vectInsTexts[pair.seqN]==1 iff the i-th text has been recovered
                if (foundSymbol == TERMINATE_CHAR)
                    vectInsTexts[pair.seqN]=1;
        
                    
                //Set newPair
                    
                //.posN
                //counters[alpha[(unsigned int)foundSymbol]]= # foundSymbol occ. in currentPile
                newPair.posN = counters[alpha[(unsigned int)foundSymbol]];

                for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {
                    //Sum #foundSymbol occ. in each pile g= 0...(currentPile-1)
                    newPair.posN = newPair.posN + tableOcc[g][alpha[(unsigned int)foundSymbol]];
                }
                //.seqN
                newPair.seqN=pair.seqN;

                #if DEBUG==1 && OMP==0
                    std::cerr << "RecoverNsymbolsReverseByVector -Result: vectInsTexts[vectTriple[k].seqN]==0 j: P[q]=" << (dataTypeNChar)newPair.posN <<  " N[q]=" << (dataTypeNSeq)newPair.seqN << std::endl << std::endl;
                #endif
                    
                newCurrentPile=alpha[(unsigned int)foundSymbol];
                
                #if OMP
                    vectPairParallelTmp[t_id][newCurrentPile].push_back (newPair);
                #else
                    vectPairTmp[newCurrentPile].push_back (newPair);
                #endif
                    
            } //end-if (vectInsTexts[vectTriple[k].seqN]==0)
            else {
                //INSERT $ in newSymb[i]
                    
                foundSymbol = TERMINATE_CHAR_LEN;

                #if USE_QS
                    newSymb[2*pair.seqN] = foundSymbol;
                    newSymb[2*pair.seqN+1] = '{';
                #else
                    newSymb[pair.seqN] = foundSymbol;
                #endif
                
                #if DEBUG==1 && OMP==0 && USE_QS==0
                    std::cerr << "RecoverNsymbolsReverseByVector +Result: vectInsTexts[" << pair.seqN << "]!=0 k=" << posInPile << " - : foundSymbol=" << (unsigned int)foundSymbol << " newSymb[pair.seqN]=" << newSymb[pair.seqN] <<  std::endl << std::endl;
                #endif
            }
                
            posInPile++;
        }//end-while ( !vectPairParallel[t_id][currentPile].empty() )
            
        fclose(InFileBWT);
        #if USE_QS
            fclose(InFileQS);
        #endif
                    
        currentPile++;

    }//end-while ( currentPile < sizeAlpha )
    
    delete [] counters;
    delete [] bufferBlock;
            
    #if DEBUG==1 && OMP==0 && USE_QS==0
        std::cerr << "NewSymbols " ;
        for (dataTypeNSeq g = 0 ; g < nText; g++) {
          std::cerr << (char)newSymb[g] << " ";
        }
    #endif
        
    #if OMP
        vectPairParallelTmp[t_id].swap(vectPairParallel[t_id]);
    #else
        vectPairTmp.swap(vectPair);
    #endif
        
    return 1;
}

//const char * fileInputPrefix = "cyc.\0";
int BCRdecode::convertFromCycFileToFastaOrFastq( string fileOutDecode)
{
    //Set maximum number of threads
    struct rlimit limit;

    /* Get max number of files in limit.rlim_cur */
    assert (getrlimit(RLIMIT_NOFILE, &limit) == 0);
    //cout << "Number of files that can be opened by this process: " << limit.rlim_cur << endl;
    
    //Set the maximum number of threads
    int maxThreads = limit.rlim_cur/(lengthRead + 1);
    if (numthreads > maxThreads){
        std::cerr << "Number of threads for convertFromCycFileToFastaOrFastq reduced to: " << maxThreads << endl;
    }
    int t_id = 0, t;
    
    //Open all cyc files
    vector < vector <FILE *> > inFilesCyc;
    
    #if OMP
        inFilesCyc.resize(maxThreads);
        for(; t_id<maxThreads; t_id++)
    #else
        inFilesCyc.resize(1);
        string fileOutName = fileOutDecode;
    #endif
        { inFilesCyc[t_id].resize(lengthRead); }
    
    #if OMP
    dataTypeNChar ind = (dataTypeNChar)ceil((long double)nText/numthreads);
    #pragma omp parallel for default(shared) firstprivate(t_id) private(t) num_threads(maxThreads)
    for (t = 0; t < numthreads; t++)
    {
        t_id = omp_get_thread_num();//id_thread
        dataTypeNSeq s_ind = t * ind; //start index
        dataTypeNChar f_ind = (t+1)*ind;  //final index
        
        if(t==numthreads-1){//last iteration
            f_ind = nText; //final index = nText for last thread
        }
        double start=omp_get_wtime();
        
        string fileOutName = fileOutDecode;
        
        if(t!=0)
            fileOutName = fileOutName + "_" + to_string(t);
    #endif
        
        //Open cyc.t_ files
        for ( dataTypelenSeq i = 0 ; i < lengthRead; i++ )
        {
            #if OMP
                string filenameIn = fileOutCyc + to_string(t) +"_" + to_string((int)(i))+".txt";
            #else
                string filenameIn = fileOutCyc + to_string((int)(i))+".txt";
            #endif
            inFilesCyc[t_id][i]= fopen(filenameIn.c_str(), "rb");
            if (inFilesCyc[t_id][i]==NULL) {
                #if OMP
                #pragma omp critical
                #endif
                {   std::cerr << "convertFromCycFileToFastaOrFastq: error opening " << filenameIn << std::endl;
                        exit (EXIT_FAILURE);
                }
            }
            fseek(inFilesCyc[t_id][i], 0, SEEK_SET);
                
        }//end-for
            
        //Open a outFile for each iteration
        std::ofstream outFile ( fileOutName );
        
        if ( outFile.is_open() == false )
        {
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "Error opening \"" << fileOutDecode << "\" file" << std::endl;
                exit ( 1 );
            }
        }

        //I must read a char for each sequence. The chars at the position i corresponds to the chars of the sequence i.
        uchar symbol, symbolQS;
        string sequence = "";
        string sequenceQS = "";
        
        #if OMP
            for ( dataTypeNSeq j = s_ind; j < f_ind; j++ )
        #else
            for ( dataTypeNSeq j = 0; j < nText ; j++ )
        #endif
        {

            #if USE_QS
                outFile << "@Seq."  << j << std::endl;
            #else
                outFile << "> Seq. "  << j << std::endl;
            #endif
            
            for ( dataTypelenSeq i = 0; i < lengthRead; i++ )
            {
                int num = fread ( &symbol, sizeof( uchar ), 1, inFilesCyc[t_id][i] );
                assert( num == 1 );
                
                #if DEBUG && OMP
                #pragma omp critical
                    {
                        cout<< "t = " << t << ", i = " << i << ", symbol = " << (char)symbol << endl;
                    }
                #endif
                
                #if USE_QS == 1
                    //assert( fread ( &symbolQS, sizeof( char ), 1, inFilesCycQual[i] ) == 1 );
                    assert( fread ( &symbolQS, sizeof( char ), 1, inFilesCyc[t_id][i] ) == 1 );
                #endif
                
                if ( (symbol!=TERMINATE_CHAR) && (symbol!=TERMINATE_CHAR_LEN) ) {
                    sequence.append ( 1, symbol );
                    #if USE_QS == 1
                        sequenceQS.append ( 1, symbolQS );
                    #endif
                }
            }
            outFile << sequence << std::endl;
            
            #if DEBUG==1
            #pragma omp critical
            {
                cerr << "Sequence: " << sequence << std::endl;
            }
            #endif
            
            sequence.clear();

            #if USE_QS == 1
                outFile << "+" << std::endl;
                outFile << sequenceQS << std::endl;
                sequenceQS.clear();
            #endif
            
        }  //end-for

        outFile.close();
        
        //Close cyc.t_ files
        for ( dataTypelenSeq i = 0; i < lengthRead; i++ )
            fclose( inFilesCyc[t_id][i] );
        
    #if OMP
        
        #pragma omp critical
        {
            std::cerr << "convertFromCycFileToFastaOrFastq: TIME THREAD " << t_id << " = " << omp_get_wtime()- start << "\n";
        }
        
    }//end-pragma
    #endif
    
    return 1;
}


FILE * BCRdecode::openFilePartialIn(dataTypedimAlpha currentPile) {
    static FILE *inFile;
    string filenameIn = "bwt_" + to_string((int)(currentPile))+ext;
    inFile = fopen(filenameIn.c_str(), "rb");
    if (inFile==NULL) {
        #if OMP
        #pragma omp critical
        #endif
        {
            std::cerr << "openFilePartialIn: file currentPile=" << (unsigned int)currentPile << ": Error opening: " << filenameIn << std::endl;
            exit (EXIT_FAILURE);
        }
    }
    return inFile;
}

dataTypeNChar BCRdecode::readOnFilePartial(uchar *buffer, dataTypeNChar toRead, FILE * InFileBWT) {
   dataTypeNChar numchar;

   numchar = fread(buffer,sizeof(uchar),toRead,InFileBWT);

   return numchar;
}

dataTypeNChar BCRdecode::writeOnFilePartial(uchar *buffer, dataTypeNChar numchar, FILE * OutFileBWT) {
   dataTypeNChar numcharWrite;
   numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
   return numcharWrite;
}

int BCRdecode::closeFilePartial(FILE * pFile) {
    fclose(pFile);
    return 1;
}

BCRdecode::~BCRdecode()
{

}
