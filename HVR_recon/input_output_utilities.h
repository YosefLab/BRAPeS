
#ifndef IOH
#define IOH

#include "declarations.h"

FragmentHashTable getFragmentHashTable(const char *fragmentFile);

std::string getVGene(const CharString& brapesIsotypeId);

/* never tested */
std::string getJGene(const CharString& brapesIsotypeId);

CharString parseVGene(CharString isotypeId);

void iterateIsotypeGenes(CharString isotypeFileName, 
						 FragmentHashTable& table);

void appendReadsAndComplements(TStringSet& readSeqs,
						  	   TStringSet& readComplementSeqs,
						       TNameSet& readIds,
						       SeqFileIn& readFile,
						       const CharString& mateId);

void loadReads(CharString readFileNameL,
		       CharString readFileNameR,
		       TStringSet& readSeqs, 
		       TNameSet& readIds);

bool parseVGenes(const char *vGene,
				 const char *filePath,
				 CharString& FR1,
				 CharString& FR2,
				 CharString& FR3,
				 CharString& CDR1,
				 CharString& CDR2);

/* not testsed */
bool parseJGenes(const char *jGene,
				 const char *filePath,
				 CharString& JRegion);



#endif






