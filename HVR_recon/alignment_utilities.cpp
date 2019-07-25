#include "alignment_utilities.h"

void getReconstructedIsotype(TString& isotype,
							 TString genomicHVR, 
							 TString reconstructedHVR) {
	TAlign align;
	resize(rows(align), 2);
	assignSource(row(align, 0), isotype);
	assignSource(row(align, 1), genomicHVR);

	unsigned int score = localAlignment(align, 
										Score<int>(1, -1, -40),
										DynamicGaps());

	for (int i = clippedBeginPosition(row(align, 0)), j = 0; j < length(reconstructedHVR); i++, j++) {
		isotype[i] = reconstructedHVR[j];
	}
}

AlignmentClass runHgAlignment(TAlign& HgAlign,
							  const double minHOverlapLength,
							  const double maxHMutationRate) {

	unsigned int HgLength = length(source(row(HgAlign, 0)));
	unsigned int readLength = length(source(row(HgAlign, 1)));

	double alignmentScore = globalAlignment(HgAlign,
										    Score<int>(1, -1, -40),
										    AlignConfig<true, true, true, true>(),
										    minHOverlapLength - readLength,
								 	    	HgLength - minHOverlapLength);

	AlignmentClass alignmentClass = getAlignmentClass(HgAlign);
	unsigned int overlapLength = getOverlapLength(HgAlign, alignmentClass);

	if (passesMutationRateTest(maxHMutationRate, alignmentScore, overlapLength)) {
		return alignmentClass;
	} else {
		return BadAlignment;
	}
}

bool passesMutationRateTest(const double mutationRate,
							const double alignmentScore,
							const unsigned int overlapLength) {
	return ((0.5 - (alignmentScore / (2 * static_cast<double>(overlapLength)))) <= mutationRate);
}

AlignmentClass getAlignmentClass(const TAlign& align) {

	unsigned int HLeadingGaps = countLeadingGaps(row(align, 0));
	unsigned int HTrailingGaps = countTrailingGaps(row(align, 0));
	unsigned int VLeadingGaps = countLeadingGaps(row(align, 1));
	unsigned int VTrailingGaps = countTrailingGaps(row(align, 1));

	if (HLeadingGaps == 0) {
		if (HTrailingGaps > 0 && VLeadingGaps > 0) {
			return HgFbOverlapAlignment;
		} else if (HTrailingGaps == 0) {
			return HgOverlapAlignment;
		} else {
			return BadAlignment;
		}
	} else if (HLeadingGaps > 0) {
		if (VTrailingGaps > 0 && HTrailingGaps == 0) {
			return FaHgOverlapAlignment;
		} else if (HTrailingGaps > 0) {
			return FaFbOverlapAlignment;
		} else {
			return BadAlignment;
		}
	} else {
		return BadAlignment;
	}
}

unsigned int getOverlapLength(const TAlign& align, AlignmentClass alignmentClass) {
	unsigned int readLength = length(source(row(align, 1)));
	if (alignmentClass == FaHgOverlapAlignment) {
		return readLength - countLeadingGaps(row(align, 0));
	} else if (alignmentClass == HgFbOverlapAlignment) {
		return readLength - countTrailingGaps(row(align, 0));
	} else if (alignmentClass == HgOverlapAlignment) {
		return readLength;
	} else if (alignmentClass == FaFbOverlapAlignment) {
		return length(source(row(align, 0)));
	} else {
		return 0;
	}
}
void alignFrameworkSeqs(const TStringSet& readSeqs,
						const TNameSet& readIds,
						const TString& Ha,
						const TString& Fg,
						const TString& Hb,
						const double minHOverlapLength,
						const double minFOverlapLength,
						const double maxHMutationRate,
						const double maxFMutationRate,
						int maximumReadsAligned,
						TStringSet& reconstructions,
						const bool isFR1,
						const bool isFR4) {
	
	if (maximumReadsAligned == -1) {
		maximumReadsAligned = (int) length(readSeqs);
	}

	RelativeReadCoordinateHashTable initialCoordinateTable;
	RelativeReadCoordinateHashTable finalCoordinateTable;
	FragmentHashTable putativeReadTable;
	int leftCoordinateBound = 0;
	int rightCoordinateBound = (int) length(Fg);

	for (unsigned int i = 0; i < maximumReadsAligned && i < length(readSeqs); i++) {
		CharString readId = readIds[i];
		TString readSeq = readSeqs[i];
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align, 0), Fg);
		assignSource(row(align, 1), readSeq);
		AlignmentClass alignmentClass = runHgAlignment(align, minFOverlapLength, maxFMutationRate); /* really is a run FgAlignment */
		
		if (alignmentClass == FaHgOverlapAlignment && Ha != "None") {
			if (isPutativeFaHgFragment(getPutativeFaHgOverlapFragment(align),
							 		   Ha, minHOverlapLength, maxHMutationRate)) {

				int initialCoordinate = (-1 * countLeadingGaps(row(align, 0)));
				initialCoordinateTable[readId] = initialCoordinate;
				if (initialCoordinate < leftCoordinateBound) {
					leftCoordinateBound = initialCoordinate;
				}
				// std::cout << "Fa -- Hg -- OVERLAP" << std::endl;
				// std::cout << align << std::endl;
				putativeReadTable[readId] = readSeq;
			}
		} else if (alignmentClass == FaHgOverlapAlignment && Ha == "None" && isFR1) {
			int initialCoordinate = (-1 * countLeadingGaps(row(align, 0)));
			initialCoordinateTable[readId] = initialCoordinate;
			if (initialCoordinate < leftCoordinateBound) {
				leftCoordinateBound = initialCoordinate;
			}
			putativeReadTable[readId] = readSeq;
		}

		else if (alignmentClass == HgFbOverlapAlignment && Hb != "None") { 
			if (isPutativeHgFbFragment(getPutativeHgFbOverlapFragment(align),
							 		   Hb, minHOverlapLength, maxHMutationRate)) {
				int initialCoordinate = length(Fg) + countTrailingGaps(row(align, 0)) - length(readSeq);
				initialCoordinateTable[readId] = initialCoordinate;
				putativeReadTable[readId] = readSeq;
				// std::cout << "Hg -- Fb -- OVERLAP" << std::endl;
				// std::cout << align << std::endl;
				if (initialCoordinate + (int)length(readSeq) > rightCoordinateBound) {
					rightCoordinateBound = initialCoordinate + (int)length(readSeq);
				}	
			} 
		} else if (alignmentClass == HgFbOverlapAlignment && Hb == "None" && isFR4) { /* test */
			//std::cout << align << std::endl;
			int initialCoordinate = length(Fg) + countTrailingGaps(row(align, 0)) - length(readSeq);
			initialCoordinateTable[readId] = initialCoordinate;
			putativeReadTable[readId] = readSeq;
			// std::cout << "Hg -- Fb -- OVERLAP" << std::endl;
			// std::cout << align << std::endl;
			if (initialCoordinate + (int)length(readSeq) > rightCoordinateBound) {
				rightCoordinateBound = initialCoordinate + (int)length(readSeq);
			}

		} else if (alignmentClass == HgOverlapAlignment) {
			int initialCoordinate;
			if (isPutativeInfixFragment(readSeq, Fg, minFOverlapLength, maxFMutationRate, initialCoordinate)) {
				//int initialCoordinate = (-1 * countLeadingGaps(row(align, 0)));
				initialCoordinateTable[readId] = initialCoordinate;
				putativeReadTable[readId] = readSeq;
			}
		}	
	}

	rightCoordinateBound = rightCoordinateBound - leftCoordinateBound;

	TStore store;
	resize(store.contigStore, 1);
	appendValue(store.contigNameStore, "ref");

	int i = 0;

	for (auto itr = initialCoordinateTable.begin(); itr != initialCoordinateTable.end(); itr++, i++) {
		int finalCoordinate = itr->second - leftCoordinateBound;
		finalCoordinateTable[itr->first] = finalCoordinate; 
		TString readSeq = putativeReadTable[itr->first];
		appendRead(store, readSeq);
		appendAlignedRead(store, i, 0, finalCoordinate, finalCoordinate + (int)length(readSeq));
	}


	/* added this!!!!! */
	// if (length(store.readSeqStore) > 1) {
	// 	ConsensusAlignmentOptions options;
 //    	options.useContigID = false;
 //    	consensusAlignment(store, options);
	// }
	/*******************/

	TString consensus = getConsensusHcRegion(putativeReadTable, finalCoordinateTable, leftCoordinateBound, rightCoordinateBound, Ha, Fg, Hb);

    //std::cout << "consensus: " << consensus << "\n"; // 

    TString Fc = infix(consensus, (-1 * leftCoordinateBound), (int)length(Fg) - leftCoordinateBound);

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), Fg);
    assignSource(row(align, 1), Fc);

    double finalScore = globalAlignment(align, Score<int>(1, -1, -40), AlignConfig<false, false, false, false>());
    if (passesMutationRateTest(maxFMutationRate, finalScore, (unsigned int)length(Fc))) {
    	std::cout << "Reconstruction Succeeded!" << std::endl;
    	appendValue(reconstructions, Fc);
    } else {
    	std::cout << "Reconstruction Failed! Using genomic sequence instead...\n";
    	appendValue(reconstructions, Fg);
    }

    //AlignedReadLayout layout;
    //layoutAlignment(layout, store);
    //printAlignment(std::cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ (int) length(consensus), 0, (int) putativeReadTable.size());	
}

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
			   TStringSet& reconstructions) {
	
	if (maximumReadsAligned == -1) {
		maximumReadsAligned = (int) length(readSeqs);
	}

	RelativeReadCoordinateHashTable initialCoordinateTable;
	RelativeReadCoordinateHashTable finalCoordinateTable;
	FragmentHashTable putativeReadTable;
	int leftCoordinateBound = 0;
	int rightCoordinateBound = (int) length(Hg);

	for (unsigned int i = 0; i < maximumReadsAligned && i < length(readSeqs); i++) {
		CharString readId = readIds[i];
		TString readSeq = readSeqs[i];
		TAlign align;
		resize(rows(align), 2);
		assignSource(row(align, 0), Hg);
		assignSource(row(align, 1), readSeq);
		AlignmentClass alignmentClass = runHgAlignment(align, minHOverlapLength, maxHMutationRate);
		if (alignmentClass == FaHgOverlapAlignment) {
			if (isPutativeFaHgFragment(getPutativeFaHgOverlapFragment(align),
							 		   Fa, minFOverlapLength, maxFMutationRate)) {

				int initialCoordinate = (-1 * countLeadingGaps(row(align, 0)));
				initialCoordinateTable[readId] = initialCoordinate;
				//std::cout << "Initial Coordinate: " << initialCoordinate << std::endl;
				if (initialCoordinate < leftCoordinateBound) {
					leftCoordinateBound = initialCoordinate;
				}
				//std::cout << "Left Coordinate Bound: " << leftCoordinateBound << std::endl;
				//std::cout << align << std::endl;
				putativeReadTable[readId] = readSeq;
			}
		} else if (alignmentClass == HgFbOverlapAlignment) {
			if (isPutativeHgFbFragment(getPutativeHgFbOverlapFragment(align),
							 		   Fb, minFOverlapLength, maxFMutationRate)) {
				int initialCoordinate = length(Hg) + countTrailingGaps(row(align, 0)) - length(readSeq);
				initialCoordinateTable[readId] = initialCoordinate;
				putativeReadTable[readId] = readSeq;
				//std::cout << "Initial Coordinate: " << initialCoordinate << std::endl;
				if (initialCoordinate + (int)length(readSeq) > rightCoordinateBound) {
					rightCoordinateBound = initialCoordinate + (int)length(readSeq);
				}
				//std::cout << "Right Coordinate Bound: " << rightCoordinateBound << std::endl;
				//std::cout << align << std::endl;	
			} 
		} else if (alignmentClass == HgOverlapAlignment) {
			int initialCoordinate;
			if (isPutativeInfixFragment(readSeq, Hg, minHOverlapLength, maxHMutationRate, initialCoordinate)) {
				//int initialCoordinate = countLeadingGaps(row(align, 0));
				initialCoordinateTable[readId] = initialCoordinate;
				putativeReadTable[readId] = readSeq;
			}
		} /* ImplementFaFbOverlapAlignment!! */
	}

	rightCoordinateBound = rightCoordinateBound - leftCoordinateBound;

	TStore store;
	resize(store.contigStore, 1);
	appendValue(store.contigNameStore, "ref");

	int i = 0;

	for (auto itr = initialCoordinateTable.begin(); itr != initialCoordinateTable.end(); itr++, i++) {
		int finalCoordinate = itr->second - leftCoordinateBound;
		finalCoordinateTable[itr->first] = finalCoordinate; 
		TString readSeq = putativeReadTable[itr->first];
		appendRead(store, readSeq);
		appendAlignedRead(store, i, 0, finalCoordinate, finalCoordinate + (int)length(readSeq));
	}

	//std::cout << "HERE1" << std::endl;
	//std::cout << "LEngth: " << length(store.readSeqStore) << std::endl;
	if (length(store.readSeqStore) > 1) {
		ConsensusAlignmentOptions options;
    	options.useContigID = false;
    	consensusAlignment(store, options);
	}
	//std::cout << "HERE2" << std::endl;

	TString consensus = getConsensusHcRegion(putativeReadTable, finalCoordinateTable, leftCoordinateBound, rightCoordinateBound, Fa, Hg, Fb);

    //std::cout << "consensus: " << consensus << "\n"; // 

    TString Hc = infix(consensus, (-1 * leftCoordinateBound), (int)length(Hg) - leftCoordinateBound);

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), Hg);
    assignSource(row(align, 1), Hc);

    double finalScore = globalAlignment(align, Score<int>(1, -1, -40), AlignConfig<false, false, false, false>());
    if (passesMutationRateTest(maxHMutationRate, finalScore, (unsigned int)length(Hc))) {
    	std::cout << "Reconstruction Succeeded!" << std::endl;
    	appendValue(reconstructions, Hc);
    } else {
    	std::cout << "Reconstruction Failed! Using genomic sequence instead...\n";
    	appendValue(reconstructions, Hg);
    }

	//AlignedReadLayout layout;
    //layoutAlignment(layout, store);
    //printAlignment(std::cout, layout, store, /*contigID=*/ 0, /*beginPos=*/ 0, /*endPos=*/ (int) length(consensus), 0, (int) putativeReadTable.size());
}


TString getConsensusHcRegion(FragmentHashTable& readTable,
							 RelativeReadCoordinateHashTable& coordinateTable,
							 const unsigned int leftCoordinateBound,
							 const unsigned int rightCoordinateBound,
							 const TString& Fa,
							 const TString& Hg,
							 const TString& Fb) {

    String<ProfileChar<Dna5>> profile;
    resize(profile, rightCoordinateBound);

    for (auto itr = readTable.begin(); itr != readTable.end(); itr++) {
    	TString readSeq = itr->second;
    	int coordinate = coordinateTable[itr->first];
    	//std::cout << "coordinate: " << coordinate << std::endl; //
    	for (int j = 0; j < length(readSeq); j++) {
    		auto ordVal = ordValue(readSeq[j]);
    		//std::cout << readSeq[j] << "\t" << ordVal << std::endl; //
    		if (ordVal != GAP_ORD_VAL) {
    			profile[coordinate + j].count[ordVal] += 1;
    			//std::cout << "profile[" << (coordinate + j) << "] = " << profile[coordinate + j] << std::endl; // 
    		}
    	}
    }

    int j = 0;
    //std::cout << "Profile length: " << length(profile) << std::endl;
    for (int i = (-1 * leftCoordinateBound); i < rightCoordinateBound; i++, j++) {
    	//std::cout << "totalCount(profile[" << i << "])" << totalCount(profile[i]) << std::endl;
    	if (totalCount(profile[i]) == 0) {
    		auto ordVal = ordValue(Hg[j]);
    		profile[i].count[ordVal] = 1;
    	}
    }

    TString consensus;

    for (int j = 0; j < length(profile); j++) {
    	appendValue(consensus, Dna5(getMaxIndex(profile[j])));
    }
    return consensus;
}

TString getPutativeFaHgOverlapFragment(const TAlign& align) {
	TString readSeq = source(row(align, 1));
	return prefix(readSeq, countLeadingGaps(row(align, 0)));

}

TString getPutativeHgFbOverlapFragment(const TAlign& align) {
	TString readSeq = source(row(align, 1));
	return suffix(readSeq, length(readSeq) - countTrailingGaps(row(align, 0)));
}

bool isPutativeInfixFragment(const TString& OverlapFragment,
						  	 const TString& Ref,
						  	 const double minRefOverlapLength,
						  	 const double maxRefMutationRate,
						  	 int& initialCoordinate) {

	unsigned int overlapLength = length(OverlapFragment);
	if (overlapLength < minRefOverlapLength || overlapLength > length(Ref)) return false;

	TAlign align;
	resize(rows(align), 2);
	assignSource(row(align, 0), Ref);
	assignSource(row(align, 1), OverlapFragment);

	double score = globalAlignment(align, Score<int>(1, -1, -40), AlignConfig<true, false, false, true>());
	initialCoordinate = countLeadingGaps(row(align, 1));

	if(passesMutationRateTest(maxRefMutationRate, score, overlapLength)) {
		//std::cout << "Putative Infix Test Alignment" << std::endl;
		//std::cout << "Score = " << score << std::endl;
		//std::cout << align << std::endl;
		return true;
	} else {
		return false;
	}
}

bool isPutativeFaHgFragment(const TString& FaHgOverlapFragment,
							const TString& Fa,
							const unsigned int minFOverlapLength,
							const double maxFMutationRate) {
	
	unsigned int overlapLength = length(FaHgOverlapFragment);
	if (overlapLength < minFOverlapLength) return false;
	
	TString pertinentFaFragment = suffix(Fa, length(Fa) - overlapLength);

	TAlign align;
	resize(rows(align), 2);
	assignSource(row(align, 0), pertinentFaFragment);
	assignSource(row(align, 1), FaHgOverlapFragment);

	double score = globalAlignment(align, Score<int>(1, -1, -40), AlignConfig<false, false, false, false>());
	unsigned int alignmentOverlapLength = getOverlapLength(align, FaHgOverlapAlignment);

	return passesMutationRateTest(maxFMutationRate, score, alignmentOverlapLength);

}

bool isPutativeHgFbFragment(const TString& HgFbOverlapFragment,
							const TString& Fb,
							const unsigned int minFOverlapLength,
							const double maxFMutationRate) {
	
	unsigned int overlapLength = length(HgFbOverlapFragment);
	if (overlapLength < minFOverlapLength) return false;
	
	TString pertinentFbFragment = prefix(Fb, overlapLength);

	TAlign align;
	resize(rows(align), 2);
	assignSource(row(align, 0), pertinentFbFragment);
	assignSource(row(align, 1), HgFbOverlapFragment);

	double score = globalAlignment(align, Score<int>(1, -1, -40), AlignConfig<false, false, false, false>());
	unsigned int alignmentOverlapLength = getOverlapLength(align, HgFbOverlapAlignment);

	if (passesMutationRateTest(maxFMutationRate, score, alignmentOverlapLength)) {
		return true;
	}
	return false;

	//return passesMutationRateTest(maxFMutationRate, score, alignmentOverlapLength);

}


/* never tested */
TString getFR4ReferenceSeq(const TString& brapesBestIsoSeq,
						   const TString& genomicJGene,
						   const unsigned int CDR3extension) {
	TAlign align;
	resize(rows(align), 2);
	assignSource(row(align, 0), brapesBestIsoSeq);
	assignSource(row(align, 1), genomicJGene);

	double score = globalAlignment(align, Score<int>(1, -1, -40), AlignConfig<true, true, true, true>());
	//TRow row = row(align, 0);
	unsigned int start = countLeadingGaps(row(align, 1));
	unsigned int end = length(brapesBestIsoSeq) - countTrailingGaps(row(align, 1));
	TString flankingCDR3 = infix(brapesBestIsoSeq, start - CDR3extension, start);
	// std::cout << genomicJGene << std::endl;
	// std::cout << infix(brapesBestIsoSeq, start, end) << std::endl;
	TString refSeq = flankingCDR3; refSeq += genomicJGene;
	return refSeq;
}











