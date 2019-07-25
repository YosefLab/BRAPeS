#include <iostream>
#include "declarations.h"
#include "input_output_utilities.h"
#include "alignment_utilities.h"
#include "parser.h"


int main(int argc, char** argv) {
  HypervariableReconstructionArguments arguments;
  HypervariableReconstructionOptions options;

  ArgumentParser::ParseResult res = parseCommandLine(arguments, options, argc, argv);

  if (res != ArgumentParser::PARSE_OK) {
    return res == ArgumentParser::PARSE_ERROR;
  }

  TStringSet readSeqs;
  TNameSet readIds;

  std::cout << "\nLoading reads...\n";
  loadReads(arguments.readFileNameL, arguments.readFileNameR,
  readSeqs, readIds);

  TStringSet reconstructions;
  TStringSet isotypeReconstructions;

  TNameSet reconstructionIds;
  TNameSet outputReconstructionIds;
  TNameSet isotypeReconstructionIds;

  SeqFileIn isos(toCString(arguments.bestIsosFileName));

  TString brapesIsotypeSeq;
  CharString brapesIsotypeId;

  CharString FR1, FR2, FR3, CDR1, CDR2, JRegion;
  int z = 0;
  while (!atEnd(isos)) {
    readRecord(brapesIsotypeId, brapesIsotypeSeq, isos);
    CharString FR1BrapesIsotypeId(brapesIsotypeId); FR1BrapesIsotypeId += " FR1";
    CharString CDR1BrapesIsotypeId(brapesIsotypeId); CDR1BrapesIsotypeId += " CDR1";
    CharString FR2BrapesIsotypeId(brapesIsotypeId); FR2BrapesIsotypeId += " FR2";
    CharString CDR2BrapesIsotypeId(brapesIsotypeId); CDR2BrapesIsotypeId += " CDR2";
    CharString FR3BrapesIsotypeId(brapesIsotypeId); FR3BrapesIsotypeId += " FR3";
    CharString FR4BrapesIsotypeId(brapesIsotypeId); FR4BrapesIsotypeId += " FR4";
    std::string vGene = getVGene(brapesIsotypeId);
    std::string jGene = getJGene(brapesIsotypeId);
    bool reconstructFR4 = false;
    std::cout << "\nWorking on " << vGene << "...\n" << std::endl;
    if (arguments.organism == "human") {
      if (parseJGenes(jGene.c_str(), "../seqs/J.human.conserved.fa", JRegion)) {
        JRegion = getFR4ReferenceSeq(brapesIsotypeSeq, JRegion, 10);
        reconstructFR4 = true;
      } else {
        //std::cout << jGene << " reference: " << JRegion << std::endl;
        //std::cout << jGene << " reference: " << JRegion << std::endl;
        std::cout << "No IMGT reference for " << jGene << ". Will skip FR4 reconstruction...\n";
      }
      if (parseVGenes(vGene.c_str(), "../seqs/human.fasta", FR1, FR2, FR3, CDR1, CDR2)) {
    
      } else {
        std::cout << "No IMGT reference for " << vGene << ". Skipping...\n";
        continue;
      } 
    } else {
      if (parseJGenes(jGene.c_str(), "../seqs/J.mouse.conserved.fa", JRegion)) {
        JRegion = getFR4ReferenceSeq(brapesIsotypeSeq, JRegion, 10);
        reconstructFR4 = true;
      } else {
        std::cout << "No IMGT reference for " << jGene << ". Will skip FR4 reconstruction...\n";
      }
      if (parseVGenes(vGene.c_str(), "../seqs/mouse.fasta", FR1, FR2, FR3, CDR1, CDR2)) {

      } else {
        std::cout << "No IMGT reference for " << vGene << ". Skipping...\n";
        continue;
      }
    }

    std::cout << "Reconstructing CDR1\n";
    alignSeqs(readSeqs,
              readIds, 
              FR1,
              CDR1,
              FR2,
              options.HMinOverlapLength,
              options.FMinOverlapLength,
              options.HMaxMutationRate,
              options.FMaxMutationRate,
              options.maximumReadsAligned,
              reconstructions);

    TString reconstructedCDR1 = reconstructions[z++];

    getReconstructedIsotype(brapesIsotypeSeq, CDR1, reconstructedCDR1);

    CharString CDR1ReconstructionId = vGene;
    CharString CDR2ReconstructionId = vGene;
    CDR1ReconstructionId += " CDR1";
    CDR2ReconstructionId += " CDR2";

    appendValue(reconstructionIds, CDR1ReconstructionId);
    appendValue(reconstructionIds, CDR2ReconstructionId);
    appendValue(outputReconstructionIds, CDR1BrapesIsotypeId);
    appendValue(outputReconstructionIds, CDR2BrapesIsotypeId);

    std::cout << "Reconstructing CDR2\n";
    alignSeqs(readSeqs,
              readIds, 
              FR2,
              CDR2,
              FR3,
              options.HMinOverlapLength,
              options.FMinOverlapLength,
              options.HMaxMutationRate,
              options.FMaxMutationRate,
              options.maximumReadsAligned,
              reconstructions);

    TString reconstructedCDR2 = reconstructions[z++];
    
    getReconstructedIsotype(brapesIsotypeSeq, CDR2, reconstructedCDR2);



    std::cout << "Reconstructing FR1\n";
    alignFrameworkSeqs(readSeqs,
                       readIds,
                       "None",
                       FR1,
                       reconstructedCDR2, /* use reconstructed version !! */
                       options.HMinOverlapLength,
                       10,
                       options.HMaxMutationRate,
                       options.FMaxMutationRate,
                       options.maximumReadsAligned,
                       reconstructions,
                       true,
                       false);

    CharString FR1ReconstructionId = vGene;
    FR1ReconstructionId += " FR1";
    appendValue(reconstructionIds, FR1ReconstructionId);
    appendValue(outputReconstructionIds, FR1BrapesIsotypeId);

    TString reconstructedFR1 = reconstructions[z++];
    getReconstructedIsotype(brapesIsotypeSeq, FR1, reconstructedFR1);

    std::cout << "Reconstructing FR2\n";
    alignFrameworkSeqs(readSeqs,
                       readIds,
                       reconstructedCDR1, /* use reconstructed version !! */
                       FR2,
                       reconstructedCDR2, /* use reconstructed version !! */
                       options.HMinOverlapLength,
                       options.FMinOverlapLength,
                       options.HMaxMutationRate,
                       options.FMaxMutationRate,
                       options.maximumReadsAligned,
                       reconstructions,
                       false,
                       false);

    CharString FR2ReconstructionId = vGene;
    FR2ReconstructionId += " FR2";
    appendValue(reconstructionIds, FR2ReconstructionId);
    appendValue(outputReconstructionIds, FR2BrapesIsotypeId);

    TString reconstructedFR2 = reconstructions[z++];
    getReconstructedIsotype(brapesIsotypeSeq, FR2, reconstructedFR2);

    std::cout << "Reconstructing FR3\n";
    alignFrameworkSeqs(readSeqs,
                       readIds,
                       reconstructedCDR2, /* use reconstructed version !! */
                       FR3,
                       "None",
                       options.HMinOverlapLength,
                       options.FMinOverlapLength,
                       options.HMaxMutationRate,
                       options.FMaxMutationRate,
                       options.maximumReadsAligned,
                       reconstructions,
                       false,
                       false);

    CharString FR3ReconstructionId = vGene;
    FR3ReconstructionId += " FR3";

    TString reconstructedFR3 = reconstructions[z++];
    getReconstructedIsotype(brapesIsotypeSeq, FR3, reconstructedFR3);

    //TString brapesIsotypeSeq; brapesIsotypeSeq += FR1; brapesIsotypeSeq += CDR1; brapesIsotypeSeq += FR2; brapesIsotypeSeq += CDR2; brapesIsotypeSeq += FR3;


    appendValue(reconstructionIds, FR3ReconstructionId);
    appendValue(outputReconstructionIds, FR3BrapesIsotypeId);
    if (reconstructFR4) {
      std::cout << "Reconstructing FR4\n";
      TString flankingCDR3 = prefix(JRegion, 10);
      TString FR4 = suffix(JRegion, 10);


      alignFrameworkSeqs(readSeqs,
                         readIds,
                         flankingCDR3, /* use reconstructed version !! */
                         FR4,
                         "None",
                         options.HMinOverlapLength,
                         10,
                         0,
                         options.FMaxMutationRate,
                         options.maximumReadsAligned,
                         reconstructions,
                         false,
                         true);

      CharString FR4ReconstructionId = vGene;
      FR3ReconstructionId += " FR4";

      TString reconstructedFR4 = reconstructions[z++];
      getReconstructedIsotype(brapesIsotypeSeq, FR4, reconstructedFR4);

      //TString brapesIsotypeSeq; brapesIsotypeSeq += FR1; brapesIsotypeSeq += CDR1; brapesIsotypeSeq += FR2; brapesIsotypeSeq += CDR2; brapesIsotypeSeq += FR3;

      appendValue(reconstructionIds, FR4ReconstructionId);
      appendValue(outputReconstructionIds, FR4BrapesIsotypeId);
    }

    appendValue(isotypeReconstructions, brapesIsotypeSeq);
    appendValue(isotypeReconstructionIds, brapesIsotypeId);



  }

  CharString outputPrefix = options.outputFileName;

  CharString outputIsotypeFileName(outputPrefix);
  outputIsotypeFileName += ".BCRs.fasta";
  CharString outputCDRsFileName(outputPrefix);
  outputCDRsFileName += ".CDR1.CDR2.fasta";

  SeqFileOut isotypeFile(toCString(outputIsotypeFileName));
  for (int i = 0; i < length(isotypeReconstructionIds); i++) {
    writeRecord(isotypeFile, isotypeReconstructionIds[i], isotypeReconstructions[i], "");
  }

  SeqFileOut reconstructionFile(toCString(outputCDRsFileName));
  for (int i = 0; i < length(reconstructions); i++) {
    writeRecord(reconstructionFile, outputReconstructionIds[i], reconstructions[i], "");
  }

  std::cout << std::endl;

	return 0;
}

