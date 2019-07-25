
#ifndef ALIGNMENT_UTILITIESH
#define ALIGNMENT_UTILITIESH

#include "declarations.h"


/*
	CDR Reconstruction Pipeline
	---------------------------

	1. Obtain all read seqs using input_output_utilities.h. For each mate pair there are 4 separate seqs used:
	   mate-pair-1 and mate-pair-1 complement; mate-pair-2 and mate-pair-2 complement. These are all passed
	   to alignSeqs as a StringSet.
	
	2. alignSeqs takes genomic reference sequences from both flanking framework regions and the CDR region. It 
	   also takes inputs from the parsed command input that specify the maximum mutation rates and the minimum
	   alignment overlap lengths.
	
	3. Loop begins on each read seq:
		
		i. 
				a. Read seqs are first aligned to the genomic reference CDR region.
				
				b. Alignment passes first test of significance if the alignment overlap
				   length is greater than minHOverlapLength and the calculated mutation rate is less
		   		   than maxHMutationRate. 
		   		
		   		c. If significant, the an AlignmentClass is assigned to the alignment depending on how
		   		   exactly the read aligned to the CDR region.
		   		
		   		d. If not signficant, toss and continue on to next read seq.
		
		ii. 
				a. If AlignmentClass is FaHgOverlapAlignment, align fragment that was part of the left 
				   overhanging gap from the CDR alignment to the the right side of the Fa region. Alignment
				   is deemed significant using the same mechanism as described above for the CDR alignment.
				
				b. If AlignmentCLass is HgFbOverlapAlignment, align fragment that was part of the right
				   overhanging gap from the CDR alignment to the left side of the Fb region. Alignment is
				   deemed significant using the same mechanism as described above.
				
				c. If AlignmentClass is FaFbOverlapAlignment, do the same as ii a. and b. except for Fa
				   and Fb.
				
				d. If AlignmentClass is HgOverlapAlignment, pass.
				
				e. If Alignment is deemed insignificant, toss and continue to next read seq.
		
		iii.	
				a. Save alignment coordinates
*/

typedef enum AlignmentClass {
	
	FaHgOverlapAlignment, 			/* read aligns across the Fa-Hg boundary */
	
	
	HgFbOverlapAlignment,			/* read aligns across the Hg-Fb boundary */
	

	FaFbOverlapAlignment,			/* read aligns across the Fa-Hg and Hg-Fb boundaries.
									   Implicitly, this means the read aligns across the 
									   entire CDR region. */	
	HgOverlapAlignment,				/* read aligns only to the CDR region. This means the read
									   is smaller than the CDR region but still aligned to it */

	BadAlignment					/* read alignment is insignificant */
} AlignmentClass;

void alignSeqs(const TStringSet& readSeqs,
			   const TNameSet& readIds,
			   const TString& Fa,
			   const TString& Hg,
			   const TString& Fb,
			   const double minHOverlapLength,
			   const double minFOverlapLength,
			   const double maxHMutationRate,
			   const double maxFMutationRate,
			   int maximumReadsAligned,
			   TStringSet& reconstructions);

void alignFrameworkSeqs(const TStringSet& readSeqs,
						const TNameSet& readIds,
						const TString& Ha,
						const TString& Fg,
						const TString& Hb,
						const double minHOverlapLength,
						const double minFOverlapLength,
						const double maxHMutationRate,
						const double maxFMutationRate,
						int maxmimumReadsAligned,
						TStringSet& reconstructions,
						const bool isFR1,
						const bool isFR4);

void getReconstructedIsotype(TString& isotype, TString genomicHVR, TString reconstructedHVR);


AlignmentClass runHgAlignment(TAlign& HgAlign,
							  const double minHOverlapLength,
							  const double maxHMutationRate);

bool passesMutationRateTest(const double mutationRate,
							const double alignmentScore,
							const unsigned int overlapLength);

unsigned int getOverlapLength(const TAlign& align, AlignmentClass alignmentClass);

AlignmentClass getAlignmentClass(const TAlign& align);

TString getPutativeFaHgOverlapFragment(const TAlign& align);

TString getPutativeHgFbOverlapFragment(const TAlign& align);

bool isPutativeFaHgFragment(const TString& FaHgOverlapFragment,
							const TString& Fa,
							const unsigned int minFOverlapLength,
							const double maxFMutationRate);

bool isPutativeHgFbFragment(const TString& HgFbOverlapFragment,
							const TString& Fb,
							const unsigned int minFOverlapLength,
							const double maxFMutationRate);

bool isPutativeInfixFragment(const TString& OverlapFragment,
						  	 const TString& Ref,
						  	 const double minRefOverlapLength,
						  	 const double maxRefMutationRate,
						  	 int& initialCoordinate);

TString getConsensusHcRegion(FragmentHashTable& readTable,
							 RelativeReadCoordinateHashTable& coordinateTable,
							 const unsigned int leftCoordinateBound,
							 const unsigned int rightCoordinateBound,
							 const TString& Fa,
							 const TString& Hg,
							 const TString& Fb);

/* never tested */
TString getFR4ReferenceSeq(const TString& brapesBestIsoSeq,
						   const TString& genomicJGene,
						   const unsigned int CDR3extension);

#endif

