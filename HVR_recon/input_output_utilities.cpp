 #include "input_output_utilities.h"



std::string getVGene(const CharString& brapesIsotypeId) {
	std::string vGene;
	for (int i = 0; i < length(brapesIsotypeId); i++) {
		if (brapesIsotypeId[i] != '.') {
			vGene.push_back(brapesIsotypeId[i]);
		} else {
			break;
		}
	}
	return vGene;
}

/* never tested */
std::string getJGene(const CharString& brapesIsotypeId) {
	std::string jGene;
	int i = 0;
	for (; brapesIsotypeId[i] != '.'; i++) {}
	i++;
	for (; i < length(brapesIsotypeId); i++) {
		if (brapesIsotypeId[i] != '.') {
			jGene.push_back(brapesIsotypeId[i]);
		} else {
			break;
		}
	}
	return jGene;
}

bool parseVGenes(const char *vGene,
				 const char *filePath,
				 CharString& FR1,
				 CharString& FR2,
				 CharString& FR3,
				 CharString& CDR1,
				 CharString& CDR2) {

	std::ifstream genomicSeqFile(toCString(getAbsolutePath(filePath)));
	std::stringstream genomicSeqFileStream;
	while (genomicSeqFile >> genomicSeqFileStream.rdbuf());

	std::stringstream regexStream;
	regexStream << ">[^>]*" << vGene << "\\*0[0-9][^>]*";
	
	std::regex searchRegex(regexStream.str());
	std::smatch infoMatch;

	const std::string seqFileStr = genomicSeqFileStream.str();
	//std::cout << seqFileStr << std::endl;

	
	std::sregex_iterator itr(seqFileStr.begin(), seqFileStr.end(), searchRegex);
	std::sregex_iterator end;
	if (itr == end) return false;


	while (itr != end) {
		for (int i = 0; i < itr->size(); i++) {
			const std::string geneInfo = (*itr)[i];
			std::regex seqRegex("\n[\nacgt\\.]*");
			std::smatch seqMatch;
			std::regex_search(geneInfo, seqMatch, seqRegex);

			std::string gene = seqMatch[0];
			gene.erase(std::remove(gene.begin(), gene.end(), '\n'), gene.end());

			std::string FR1str;
      		std::string CDR1str;
      		std::string FR2str;
      		std::string CDR2str;
      		std::string FR3str;

      		int j = 0;
      		
      		for (; j < 78; j++) {
        		FR1str.push_back(gene[j]);
    		}

    		for (; j < 114; j++) {
        		CDR1str.push_back(gene[j]);
      		}

      		for (; j < 165; j++) {
        		FR2str.push_back(gene[j]);
      		}

      		for (; j < 195; j++) {
        		CDR2str.push_back(gene[j]);
      		}

      		for (; j < 309; j++) {
        		FR3str.push_back(gene[j]);
      		}

      		FR1str.erase(std::remove(FR1str.begin(), FR1str.end(), '.'), FR1str.end());
      		CDR1str.erase(std::remove(CDR1str.begin(), CDR1str.end(), '.'), CDR1str.end());
      		FR2str.erase(std::remove(FR2str.begin(), FR2str.end(), '.'), FR2str.end());
      		CDR2str.erase(std::remove(CDR2str.begin(), CDR2str.end(), '.'), CDR2str.end());
      		FR3str.erase(std::remove(FR3str.begin(), FR3str.end(), '.'), FR3str.end());

      		FR1 = FR1str;
      		CDR1 = CDR1str;
      		FR2 = FR2str;
      		CDR2 = CDR2str;
      		FR3 = FR3str;

      		toUpper(FR1);
      		toUpper(CDR1);
      		toUpper(FR2);
      		toUpper(CDR2);
      		toUpper(FR3);   /* Something weird going on here, make this part more logical */ /* I think it's been fixed */

		}
		break;
	}

	// std::cout << FR1 << '\n';
	// std::cout << FR2 << '\n';
	// std::cout << FR3 << '\n';
	// std::cout << CDR1 << '\n';
	// std::cout << CDR2 << '\n';

	return true;
}
/* not testsed */
bool parseJGenes(const char *jGene,
				 const char *filePath,
				 CharString& JRegion) {

	std::ifstream genomicSeqFile(toCString(getAbsolutePath(filePath)));
	std::stringstream genomicSeqFileStream;
	while (genomicSeqFile >> genomicSeqFileStream.rdbuf());

	std::stringstream regexStream;
	regexStream << ">[^>]*" << jGene << "\\*0[0-9][^>]*";

	std::regex searchRegex(regexStream.str());
	std::smatch infoMatch;

	const std::string seqFileStr = genomicSeqFileStream.str();
	
	std::sregex_iterator itr(seqFileStr.begin(), seqFileStr.end(), searchRegex);
	std::sregex_iterator end;
	if (itr == end) return false;
	

	while (itr != end) {
		for (int i = 0; i < itr->size(); i++) {
			
			const std::string geneInfo = (*itr)[i];
			std::regex seqRegex("\n[\nacgt\\.]*");
			std::smatch seqMatch;
			std::regex_search(geneInfo, seqMatch, seqRegex);
			
			std::string gene = seqMatch[0];
			gene.erase(std::remove(gene.begin(), gene.end(), '\n'), gene.end());

			std::string jRegion;
			

      		int j = 0;
      		
      		for (; j < gene.size(); j++) {
      			jRegion.push_back(gene[j]);
      		}

      		jRegion.erase(std::remove(jRegion.begin(), jRegion.end(), '.'), jRegion.end()); 		

      		JRegion = jRegion;		

      		toUpper(JRegion);
		}
		break;
	}

	return true;	

}



void appendReadsAndComplements(TStringSet& readSeqs,
						  	   TStringSet& readComplementSeqs,
						       TNameSet& readIds,
						       SeqFileIn& readFile,
						       const CharString& mateId) {
	CharString id;
	TString readSeq;

	while (!atEnd(readFile)) {
		readRecord(id, readSeq, readFile);
		append(id, mateId);
		appendValue(readSeqs, readSeq);
		appendValue(readIds, id);

		TString readComplement(readSeq); reverseComplement(readComplement);


		appendValue(readComplementSeqs, readComplement);
	}
}

void loadReads(CharString readFileNameL,
		       CharString readFileNameR,
		       TStringSet& readSeqs, 
		       TNameSet& readIds) {
	
	TStringSet readComplementSeqs;

	SeqFileIn readFileL(toCString(readFileNameL));
	SeqFileIn readFileR(toCString(readFileNameR));

	appendReadsAndComplements(readSeqs, readComplementSeqs, readIds, readFileL, ";Left-Mate");
	appendReadsAndComplements(readSeqs, readComplementSeqs, readIds, readFileR, ";Right-Mate");

	TNameSet readComplementIds;
	for (unsigned int i = 0; i < length(readIds); i++) {
		CharString readComplementId(readIds[i]);
		append(readComplementId, ";Complement");
		appendValue(readComplementIds, readComplementId);
	}

	append(readSeqs, readComplementSeqs);
	append(readIds, readComplementIds);
}