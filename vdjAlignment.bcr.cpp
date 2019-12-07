/*
 * vdjAlignment.cpp
 *
 *  Created on: Oct 3, 2014
 *      Author: afiks
 */

#include <iostream>
#include <string>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <vector>
#include <stdlib.h>

using namespace std;
using namespace seqan;

typedef String<char> TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef typename Cols<TAlign>::Type TCols;
typedef typename Row<TAlign>::Type TRow;
typedef typename Iterator<TRow>::Type TRowsIterator;

TAlign * alignToSegment( Dna5String readSeq, Dna5String cSegment, TAlign * alignment, TAlign * revCompAlign, bool isV, int scoreTh) {
	AlignConfig<true, true, true, false> alignConfigJ;
	AlignConfig<false, true, true, true > alignConfigV;
	TAlign alignmet;
	TAlign revCompAlignment;
	TSequence seqRead = readSeq;
	TSequence revCompRead = readSeq;
	reverseComplement(revCompRead);
	TSequence seqTrans = cSegment;
	resize(rows(*alignment), 2);
	assignSource(row(*alignment, 0), seqRead);
	assignSource(row(*alignment, 1), seqTrans);
	Score<int> scoringScheme(1, -1, -4, -20);
	int score = 0;
	int scoreRevComp = 0;
	if (isV==true) {
		score = globalAlignment(*alignment, scoringScheme, alignConfigV);
	} else {
		score = globalAlignment(*alignment, scoringScheme, alignConfigJ);
	}
/*	// Do the same only for the reverse complement:
	resize(rows(*revCompAlign), 2);
	assignSource(row(*revCompAlign, 0), revCompRead);
	assignSource(row(*revCompAlign, 1), seqTrans);
	if (isV==true) {
		scoreRevComp = globalAlignment(*revCompAlign, scoringScheme, alignConfigV);
	} else {
		scoreRevComp = globalAlignment(*revCompAlign, scoringScheme, alignConfigJ);
	}*/
	if (score >= scoreTh) {
		//cout << score << endl;
		//cout << *alignment << endl;
		return alignment;
	} else {
		return NULL;
	}
}




Dna5String extendSegment(vector<TAlign> readsAlign, bool isV, Dna5String former, int maxLength) {
	Dna5String fString = "";
	int formLength = length(former);
	//cout << formLength << endl;
	int posMat[4][maxLength];
	for (int i = 0; i<4; ++i) {
		for (int j = 0; j < maxLength; ++j) {
			posMat[i][j] = 0;
		}
	}
	for (vector<TAlign>::iterator alignment = readsAlign.begin(); alignment != readsAlign.end(); ++alignment) {
		TRow readRow = row(*alignment, 0);
		TRow segRow = row(*alignment, 1);
		TRowsIterator readIt = begin(readRow);
		TRowsIterator segIt = begin(segRow);
		if (isV == false) {
			goEnd(readIt);
			goEnd(segIt);
		} else {
			goBegin(readIt);
			goBegin(segIt);

		}
		int pos = 0;
		Dna5String currBase;
		while (((isV == true) & (!atEnd(readIt))) | ((isV == false) & (!atBegin(readIt)))) {
			if (isV == false) {
				goPrevious(readIt);
				goPrevious(segIt);
			}
			if (getValue(readIt) == '-') {
				currBase = 'N';
			} else {
				currBase = getValue(readIt);
			}
			if (currBase == 'A') {
				++posMat[0][pos];
			} else if (currBase == 'G') {
				++posMat[1][pos];
			} else if (currBase == 'C') {
				++posMat[2][pos];
			} else if (currBase == 'T') {
				++posMat[3][pos];
			}
			if (isV == true) {
				goNext(segIt);
				goNext(readIt);
			}
			++pos;
		}
	}
	int count = 0;
	int countRev = maxLength -1;
	int posIter = 1;
	int foundAlign = 0;
	bool lastZero = false;
	int jj = 0;
	// Final is the first non-zero position of a matrix
	while ((lastZero == false) & (jj<maxLength)) {
		int maxC = posMat[0][jj];
		if (posMat[1][jj] > maxC) {
			maxC = posMat[1][jj];
		}
		if (posMat[2][jj] > maxC) {
			maxC = posMat[2][jj];
		}
		if (posMat[3][jj] > maxC) {
			maxC = posMat[3][jj];
		}
		if (maxC > 0) {
			lastZero = true;
		}
		jj++;
	}
	int final = jj - 1;
	while (((isV == true) & (count < maxLength)) | ((isV == false) & (countRev >= 0))) {
		int curr;
		if (isV == true) {
			curr = count;
		} else {
			curr = countRev;
		}
		//cout << former << endl;
		//cout << curr << endl;
		int max = posMat[0][curr];
		int maxPos = 0;
		if (posMat[1][curr] > max) {
			max = posMat[1][curr];
			maxPos = 1;
		}
		if (posMat[2][curr] > max) {
			max = posMat[2][curr];
			maxPos = 2;
		}
		if (posMat[3][curr] > max) {
			max = posMat[3][curr];
			maxPos = 3;
		}
		Dna5String toAdd;
		//cout << max << endl;
		//cout << "gen value" << endl;
		//if (curr < formLength) {
		//	cout << getValue(former,curr) << endl;
		//}
		//cout << "end get value" << endl;
		//if (isV == true) {
		//cout << former << endl;
		//cout << curr << endl;
		//cout << formLength << endl;
		//cout << maxLength << endl;

		if (max > 0) {
			foundAlign ++;
			if (maxPos == 0) {
				toAdd = "A";
			} else if (maxPos == 1) {
				toAdd = "G";
			} else if (maxPos == 2) {
				toAdd = "C";
			} else if (maxPos == 3) {
				toAdd = "T";
			}
		} else {
			if (isV == true) {
				if (foundAlign < 1) {
				//if (posIter < 25) {
					toAdd = getValue(former,curr);
				}
			} else {
				if (curr < final) {
                //if (curr < 25) {
					//cout << "inside the if" << endl;
					//cout << curr << endl;
					toAdd = getValue(former, formLength - curr - 1);
				}
			}
		}
		//if (isV == false) {
		//	cout << former << endl;
		//	cout << max << endl;
		//	cout << posIter << endl;
		//	cout << formLength << endl;
		//	cout << curr << endl;
		//	cout << toAdd << endl;
		//}
		posIter++;
		//} else {
		//	cout << former << endl;
		//
		//	if (max > 0) {
		//		if (maxPos == 0) {
		//			toAdd = "A";
		//		} else if (maxPos == 1) {
		//			toAdd = "G";
		//		} else if (maxPos == 2) {
		//			toAdd = "C";
		//		} else if (maxPos == 3) {
		//			toAdd = "T";
		//		}
		//	}
		//	cout << toAdd << endl;
		//}
		if (isV == true) {
			count++;
		} else {
			countRev--;
		}
		//if (isV == false) {
		//	cout << toAdd << endl;
		//}
		append(fString, toAdd);
	}
	return fString;
}




int segmentsOverlap(Dna5String vSeg, Dna5String jSeg, int overlapTh) {
	int len = length(jSeg);
	int lenV = length(vSeg);
	if (len > lenV) {
		len = lenV;
	}
	for (int i = len-1; i>=overlapTh; --i) {
		if (prefix(jSeg, i) == suffix(vSeg, lenV-i)) {
			return i;
		}
	}
	return -1;

}

Dna5String expandTranscriptSequnece(Dna5String transcriptSeq, const char * readsFile, int numIterations, int scoreTh, int overlapTh) {
	bool conv = false;
	int convCount = 0;
	CharString readId;
	Dna5String readSeq;
	Dna5String finalTranscript;
	int transcriptLength = length(transcriptSeq);
	int jSegmentStart = transcriptLength/2;
	Dna5String vSegment = prefix(transcriptSeq, jSegmentStart);
	Dna5String jSegment = suffix(transcriptSeq, jSegmentStart);
	while (conv == false) {
		SequenceStream readsStream(readsFile);
		if (!isGood(readsStream)) {
			cerr << "ERROR! Could not open sequencing file\n";
			return 1;
		}

		convCount++;
		vector<TAlign> vSegmentsAlignment;
		vector<TAlign> jSegmentsAlignment;
		int maxLengthVseg = 0;
		int maxLengthJseg = 0;
		while (!atEnd(readsStream)) {
			if (readRecord(readId, readSeq, readsStream) != 0) {
				cerr << "ERROR! Could not read sequences from sequencing file\n";
				return 1;
			}
			// Perform sequence alignment for the V segment
			TAlign align;
			TAlign revCompAlign;
			TAlign * vAlignment = alignToSegment(readSeq, vSegment, &align, &revCompAlign, true, scoreTh);
			if (vAlignment != NULL) {
				//cout << *vAlignment << endl;
				vSegmentsAlignment.push_back(*vAlignment);
				int currLength = length(row(*vAlignment, 0));
				if (currLength > maxLengthVseg) {
					maxLengthVseg = currLength;
				}
			}
			// Now align to the J segments
			TAlign jAlign;
			TAlign jRevCompAlign;
			TAlign * jAlignment = alignToSegment(readSeq, jSegment, &jAlign, &jRevCompAlign, false, scoreTh);
			if (jAlignment != NULL) {
				int currLength = length(row(*jAlignment, 0));
				if (currLength > maxLengthJseg) {
					maxLengthJseg = currLength;
					}
				//cout << *jAlignment << endl;
				jSegmentsAlignment.push_back(*jAlignment);
			}
		}
		Dna5String newJseg = jSegment;
		if (jSegmentsAlignment.empty() == false) {
			newJseg = extendSegment(jSegmentsAlignment, false, jSegment, maxLengthJseg);
		}
		Dna5String newVseg = vSegment;
		if (vSegmentsAlignment.empty() == false) {
			newVseg = extendSegment(vSegmentsAlignment, true, vSegment, maxLengthVseg);
		}
		int pos = segmentsOverlap(newVseg, newJseg, overlapTh);

		if (pos != -1) {
			conv = true;
			finalTranscript = newVseg;
			finalTranscript += suffix(newJseg, pos);
		} else {
			//cout << "newVseg : " << newVseg << endl;
			//cout << "vSegment: " << vSegment << endl;
			//cout << "newJseg : " << newJseg << endl;
			//cout << "jSegment: " << jSegment << endl;
 			if ((newVseg == vSegment) && (newJseg == jSegment)) {
				finalTranscript = newVseg;
				finalTranscript += "______";
				finalTranscript += newJseg;
				conv = true;
			}
 			if (convCount > numIterations) {
 				finalTranscript = newVseg;
 				finalTranscript += "_________________________________________";
 				finalTranscript += newJseg;
 				conv = true;
 			}
			vSegment = newVseg;
			jSegment = newJseg;
		}
	}
	//cout << "Converged after " << convCount << " iterations\n" ;
	return finalTranscript;
}




/*
 *
 */
int main(int argc, char *argv[])
{
	if (length(argv) != 7 ) {
		cerr << "USAGE: vdjAlignment <readsFile> <transcriptsFile> <outputFile> <numIterations> <thresholdScore> <overlap>";
				return 1;
	}
	char * readsFile = argv[1];
	char * transcriptFile = argv[2];
	char * outputFile = argv[3];
	int numIterations = atoi(argv[4]);
	int scoreTh = atoi(argv[5]);
	int overlapTh = atoi(argv[6]);
	SequenceStream transcriptStream(transcriptFile);
	SequenceStream outputStream(outputFile, SequenceStream::WRITE);
	if (!isGood(outputStream)) {
		cerr << "ERROR! Could not create output file\n";
		return 1;
	}
	if (!isGood(transcriptStream)) {
		cerr << "ERROR! Could not open transcripts file\n";
		return 1;
	}
	CharString transcriptId;
	Dna5String transcriptSeq;
	while (!atEnd(transcriptStream)) {
		if (readRecord(transcriptId, transcriptSeq, transcriptStream) != 0) {
			cerr << "ERROR! Could not read sequences from sequencing file\n";
			return 1;
		}
		//cout << transcriptId << endl;
		//cout << transcriptSeq << endl;
		Dna5String finalTranscript = expandTranscriptSequnece(transcriptSeq, readsFile, numIterations, scoreTh, overlapTh);
		if (writeRecord(outputStream, transcriptId, finalTranscript) != 0) {
			cerr << "ERROR: Could not write transcript to output file\n";
			return 1;
		}
	}
	return 0;
}






