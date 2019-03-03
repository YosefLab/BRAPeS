#include <iostream>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/consensus.h>
#include <seqan/store.h>
#include <seqan/modifier.h>
#include <vector>



using namespace seqan;

typedef Dna5String TString;
typedef Align<TString> TAlign;

TString eliminateFramework(TString genomicHVR, int frameworkExtension) {
	TString returnHVR = "";
	for (int i = frameworkExtension; i < length(genomicHVR) - frameworkExtension; i++) {
		append(returnHVR, genomicHVR[i]);
	}
	return returnHVR;
}

TString reconstructHVRManually(std::deque<TString> reads_deque, TString genomicHVR, int frameworkExtension) {
	
	std::random_shuffle(reads_deque.begin(), reads_deque.end());
	std::queue<TString> reads(reads_deque);

	TString consensus = reads.front();
	reads.pop();

	TString genomic = eliminateFramework(genomicHVR, frameworkExtension);
	//std::cout << genomic << std::endl;

	TAlign align;
	resize(rows(align), 2);
	assignSource(row(align, 0), consensus);

	while (!reads.empty()) {
		TString seq = reads.front();
		reads.pop();

		assignSource(row(align, 1), seq);

		int score = globalAlignment(align, Score<int>(1, -1, -4, -20), AlignConfig<true, true, true, true>());
		// std::cout << align << std::endl;
	
		TString leftGap;
		TString overlap;
		TString rightGap;

		unsigned int consensusLeadingGaps = countLeadingGaps(row(align, 0));
		unsigned int consensusTrailingGaps = countTrailingGaps(row(align, 0));
		unsigned int seqLeadingGaps = countLeadingGaps(row(align, 1));
		unsigned int seqTrailingGaps = countTrailingGaps(row(align, 1));

		if (consensusLeadingGaps > 0) {
			leftGap = infix(seq, 0, consensusLeadingGaps);
			overlap = infix(consensus, 0, length(consensus) - seqTrailingGaps);
			rightGap = infix(consensus, length(consensus) - seqTrailingGaps, length(consensus));
			
			if (score == length(overlap)) {
				consensus = "";
				append(consensus, leftGap);
				append(consensus, overlap);
				append(consensus, rightGap);

				// std::cout << "Left Gap: " << leftGap << std::endl;
				// std::cout << "Overlap: " << overlap << std::endl;
				// std::cout << "Right Gap: " << rightGap << std::endl;
				assignSource(row(align, 0), consensus);

			}

		} else if (consensusTrailingGaps > 0) {
			leftGap = infix(consensus, 0, seqLeadingGaps);
			overlap = infix(seq, 0, length(seq) - consensusTrailingGaps);
			rightGap = infix(seq, length(seq) - consensusTrailingGaps, length(seq));

			if (score == length(overlap)) {
				consensus = "";
				append(consensus, leftGap);
				append(consensus, overlap);
				append(consensus, rightGap);
				assignSource(row(align, 0), consensus);
				// std::cout << "Left Gap: " << leftGap << std::endl;
				// std::cout << "Overlap: " << overlap << std::endl;
				// std::cout << "Right Gap: " << rightGap << std::endl;
			}

		} else if (consensusLeadingGaps == 0 && consensusTrailingGaps == 0) {
			//std::cout << consensus << std::endl;
		}
    }

    assignSource(row(align, 0), consensus);
    assignSource(row(align, 1), genomic);

    int score = globalAlignment(align, Score<int>(1, -1, -4, -20), AlignConfig<true, true, true, true>());

    std::cout << score << std::endl;
    std::cout << align << std::endl;

    TString HVR = infix(consensus, countLeadingGaps(row(align, 1)), length(consensus) - countTrailingGaps(row(align, 1)));
    return HVR;
}

TString repeatedVerification(std::deque<TString>& reads, TString genomicHVR, int frameworkExtension, int repeats) {
	std::vector<TString> reconstructions;
	for (int i = 0; i < repeats; i++) {
		TString reconstruction = reconstructHVRManually(reads, genomicHVR, frameworkExtension);
		reconstructions.push_back(reconstruction);
	}

	std::map<TString, int> freq_map;
	for (auto const & x : reconstructions)
    	++freq_map[x];
	
	std::vector<TString> uv;
	std::vector<int> freq_uv;
	
	for (auto const & p : freq_map)
	{
    	uv.push_back(p.first);
    	freq_uv.push_back(p.second);
	}

	int max_count = -99999;
	TString reconstructedHVR = "";

	for (int i = 0; i < uv.size(); i++) {
		if (freq_uv[i] > max_count) {
			max_count = freq_uv[i];
			reconstructedHVR = uv[i];
		}
	}

	double best_frequency = (double) max_count / reconstructions.size();
	std::cout << "Best Reconstruction Frequency: " << best_frequency << std::endl;
	//std::cout << "Best Reconstruction: " << reconstructedHVR << std::endl;
	return reconstructedHVR;
}

TString reconstructHVR(FragmentStore<>& store, TString genomicHVR, int readCount, int frameworkExtension) {

	genomicHVR = eliminateFramework(genomicHVR, frameworkExtension);

	std::vector<char> reconstruction;
	for (int i = 0; i < length(genomicHVR); i++) {
		reconstruction.push_back(genomicHVR[i]);
	}


	ConsensusAlignmentOptions options;
    options.useContigID = false;
    consensusAlignment(store, options);
    CharString HVR = store.contigStore[0].seq;
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    std::cout << "Store Length: " << length(store) << std::endl;
    int HVRLength = length(HVR);
    printAlignment(std::cout, layout, store, 0, 0, HVRLength, 0, readCount);
    
    TAlign align;
	resize(rows(align), 2);
	assignSource(row(align, 0), genomicHVR);
	assignSource(row(align, 1), HVR);
	Score<int> scoringScheme(1, -1, -4, -20);
	localAlignment(align, scoringScheme, DynamicGaps());

	TString reconstructedHVR = "";
	
	int start = clippedBeginPosition(row(align, 1));
	int end = clippedEndPosition(row(align, 1));


	for (int i = start; i < start + length(genomicHVR) && i < length(HVR); i++) {
		//append(reconstructedHVR, HVR[i]);
		reconstruction[i - start] = HVR[i];
	}

	for (int i = 0; i < reconstruction.size(); i++) {
		append(reconstructedHVR, reconstruction[i]);
	}


	// int variableLength = length(genomicHVR) - (2 * frameworkExtension);
	// int overExtension = length(HVR) - length(reconstructedHVR) - start;
	// std::cout << "Variable Length: " << variableLength << std::endl;
	// std::cout << "Over extension: " << overExtension << std::endl;
	// TString adjustedHVR = "";

	// if (overExtension > 0) {
	// 	for (int i = frameworkExtension; i < length(reconstructedHVR) - frameworkExtension; i++) {
	// 		append(adjustedHVR, reconstructedHVR[i]);
	// 		std::cout << adjustedHVR << std::endl;
	// 	}
	// } else if (overExtension <= 0) {
	// 	for (int i = frameworkExtension; i < length(reconstructedHVR); i++) {
	// 		append(adjustedHVR, reconstructedHVR[i]);
	// 		std::cout << adjustedHVR << std::endl;
	// 	}
	// } 

	return reconstructedHVR;
}

void runGlobalReadAlignment(TAlign* align, TString* read, TString genomicHVR, int& score) {
	int forwardScore;
	int reverseScore;

	TAlign forwardAlign;
	TAlign reverseAlign;

	TString forwardRead = *read;
	TString reverseRead = forwardRead;
	reverseComplement(reverseRead);

	Score<int> scoringScheme(1, -1, -4, -20);

	resize(rows(forwardAlign), 2);
	assignSource(row(forwardAlign, 0), genomicHVR);
	assignSource(row(forwardAlign, 1), forwardRead); 

	forwardScore = globalAlignment(forwardAlign, scoringScheme, AlignConfig<true, true, true, true>(), LinearGaps());

	resize(rows(reverseAlign), 2);
	assignSource(row(reverseAlign, 0), genomicHVR);
	assignSource(row(reverseAlign, 1), reverseRead);

	reverseScore = globalAlignment(reverseAlign, scoringScheme, AlignConfig<true, true, true, true>(), LinearGaps());

	if (forwardScore > reverseScore) {
		*align = forwardAlign;
		*read = forwardRead;
		score = forwardScore;
	} else {
		*align = reverseAlign;
		*read = reverseRead;
		score = reverseScore;
	}
}



void runLocalReadAlignment(TAlign* align, TString* read, TString genomicHVR, int& score) {

	int forwardScore;
	int reverseScore;

	TAlign forwardAlign;
	TAlign reverseAlign;

	TString forwardRead = *read;
	TString reverseRead = forwardRead;
	reverseComplement(reverseRead);

	Score<int> scoringScheme(1, -1, -4, -20);
	resize(rows(forwardAlign), 2);
	assignSource(row(forwardAlign, 0), genomicHVR);
	assignSource(row(forwardAlign, 1), forwardRead);

	forwardScore = localAlignment(forwardAlign, scoringScheme, DynamicGaps());

	resize(rows(reverseAlign), 2);
	assignSource(row(reverseAlign, 0), genomicHVR);
	assignSource(row(reverseAlign, 1), reverseRead);

	reverseScore = localAlignment(reverseAlign, scoringScheme, DynamicGaps());

	if (forwardScore > reverseScore) {
		*align = forwardAlign;
		*read = forwardRead;
		score = forwardScore;
	} else {
		*align = reverseAlign;
		*read = reverseRead;
		score = reverseScore;
	}
}

bool containsN(CharString read) {
	for (int i = 0; i < length(read); i++) {
		if (read[i] == 'N') {
			return true;
		}
	}
	return false;
}

void runPipe(CharString bamFileName, TString genomicCDR1, TString genomicCDR2, TString* reconstructedCDR1, TString* reconstructedCDR2, 
	int thresholdScore, int minReads, int maxReads, int readsBeforeMovingOn, CharString useAllReadsString, int frameworkExtension) {

	bool useAllReads;
	if (useAllReadsString == "All") {
		useAllReads = true;
	} else if (useAllReadsString == "Unmapped") {
		useAllReads = false;
	} else {
		throw Exception();
	}

	FragmentStore<> CDR1Store;
	FragmentStore<> CDR2Store;

	// std::queue<TString> CDR1queue;
	// std::queue<TString> CDR2queue;
	std::deque<TString> CDR1deque;
	std::deque<TString> CDR2deque;
	
	int CDR1ReadCount = 0;
	int CDR2ReadCount = 0;
	int unmappedReadsCount = 0;

	bool searchingForCDR1Reads = true;
	bool searchingForCDR2Reads = true;

	bool foundReadsForCDR1 = false;
	bool foundReadsForCDR2 = false;

	BamFileIn bamFileIn(toCString(bamFileName));
	BamHeader header;
	readHeader(header, bamFileIn);
	BamAlignmentRecord record;
	while (!atEnd(bamFileIn) && (searchingForCDR1Reads || searchingForCDR2Reads)) {
		readRecord(record, bamFileIn);
		if (useAllReads || (hasFlagUnmapped(record) && !useAllReads)) {

			TString read = record.seq;
			if(containsN(read)) {
				continue;
			}

			TAlign align;
			int score;

			if (searchingForCDR1Reads) {
				runLocalReadAlignment(&align, &read, genomicCDR1, score);
				if (score >= thresholdScore) {
					appendRead(CDR1Store, read);
					CDR1deque.push_back(read);
					std::cout << "CDR1 " << read << " " << record.qName << std::endl;
					CDR1ReadCount++;
				}
			}	
			
			if (searchingForCDR2Reads) {
				runLocalReadAlignment(&align, &read, genomicCDR2, score);
				if (score >= thresholdScore) {
					appendRead(CDR2Store, read);
					CDR2deque.push_back(read);
					std::cout << "CDR2 " << read << " " << record.qName << std::endl;
					CDR2ReadCount++;
				}
			}

			if (unmappedReadsCount > readsBeforeMovingOn) {
				searchingForCDR1Reads = false;
				searchingForCDR2Reads = false;
			} else {
				if (CDR1ReadCount > maxReads) {
					searchingForCDR1Reads = false;
				}
				if (CDR2ReadCount > maxReads) {
					searchingForCDR2Reads = false;
				}
			}
			
			foundReadsForCDR1 = (CDR1ReadCount >= minReads);
			foundReadsForCDR2 = (CDR2ReadCount >= minReads);
			unmappedReadsCount++;	
		}
	}
	if (foundReadsForCDR1) {
		std::cout << "Using " << CDR1ReadCount << " reads to reconstruct CDR1..." << std::endl;
		//*reconstructedCDR1 = reconstructHVR(CDR1Store, genomicCDR1, CDR1ReadCount, frameworkExtension);
		*reconstructedCDR1 = repeatedVerification(CDR1deque, genomicCDR1, frameworkExtension, 10);
		std::cout << "Reconstructed CDR1: " << *reconstructedCDR1 << std::endl;
	} else {
		std::cout << "Unable to reconstruct CDR1, there were only " << CDR1ReadCount << 
		" reads that aligned above the min-reads threshold (Current min-reads threshold is " << minReads << ")" << std::endl;
		std::cout << "Run again with a different score threshold or min-reads threshold to get reconstruction" << std::endl;
 	}

 	std::cout << std::endl;

	if (foundReadsForCDR2) {
		std::cout << "Using " << CDR2ReadCount << " reads to reconstruct CDR2..." << std::endl;
		//*reconstructedCDR2 = reconstructHVR(CDR2Store, genomicCDR2, CDR2ReadCount, frameworkExtension);
		*reconstructedCDR2 = repeatedVerification(CDR2deque, genomicCDR2, frameworkExtension, 10);
		std::cout << "Reconstructed CDR2: " << *reconstructedCDR2 << std::endl;
 	} else {
		std::cout << "Unable to reconstruct CDR2, there were only " << CDR2ReadCount << 
		" reads that aligned above the min-reads threshold (Current min-reads threshold is " << minReads  << ")" << std::endl;
		std::cout << "Run again with a different score threshold or min-reads threshold to get reconstruction" << std::endl;
	}
}

void correctBestIso(TString* bestIso, TString genomicCDR1, TString genomicCDR2, TString reconstructedCDR1, TString reconstructedCDR2) {
	
	TAlign align;
	resize(rows(align), 2);
	Score<int> scoringScheme(1, -1, -4, -20);
	assignSource(row(align, 0), *bestIso);

	if (reconstructedCDR1 != "None") {
		assignSource(row(align, 1), genomicCDR1);
		localAlignment(align, scoringScheme, DynamicGaps());
		for (int i = clippedBeginPosition(row(align, 0)), j = 0; j < length(reconstructedCDR1); i++, j++) {
			(*bestIso)[i] = reconstructedCDR1[j];
		}
	}

	if (reconstructedCDR2 != "None") {
		assignSource(row(align, 1), genomicCDR2);
		localAlignment(align, scoringScheme, DynamicGaps());
		for (int i = clippedBeginPosition(row(align, 0)), j = 0; j < length(reconstructedCDR2); i++, j++) {
			(*bestIso)[i] = reconstructedCDR2[j];
		}
	}
}


int main(int argc, char** argv) {

	CharString bamFileName = argv[1];
	TString genomicCDR1 = argv[2];
	TString genomicCDR2 = argv[3];
	int thresholdScore = atoi(argv[4]);
	int minReads = atoi(argv[5]);
	int maxReads = atoi(argv[6]);
	int readsBeforeMovingOn = atoi(argv[7]);
	TString bestIso = argv[8];
	CharString useAllReadsString = argv[9];
	int frameworkExtension = atoi(argv[10]);

	TString reconstructedCDR1 = "None";
	TString reconstructedCDR2 = "None";

	//std::cout << "HERE" << std::endl;

	// std::cout << "genomicCDR1: " << genomicCDR1 << std::endl;
	// std::cout << "genomicCDR2: " << genomicCDR2 << std::endl;
	
	runPipe(bamFileName, genomicCDR1, genomicCDR2, &reconstructedCDR1, &reconstructedCDR2, thresholdScore, minReads, maxReads, readsBeforeMovingOn, useAllReadsString, frameworkExtension);

	SeqFileOut HVRFileOut("hvr.reconstructions.temp.fasta");
	SeqFileOut BCRFileOut("bcr.reconstructions.temp.fasta");
	
	if (reconstructedCDR1 != "None") {
		writeRecord(HVRFileOut, "CDR1", reconstructedCDR1, "");
	}

	if (reconstructedCDR2 != "None") {
		writeRecord(HVRFileOut, "CDR2", reconstructedCDR2, "");
	}

	correctBestIso(&bestIso, genomicCDR1, genomicCDR2, reconstructedCDR1, reconstructedCDR2);
	writeRecord(BCRFileOut, "", bestIso, "");

	return 0;
}

