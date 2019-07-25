
#ifndef PARSERH
#define PARSERH

#include <iostream>
#include <seqan/arg_parse.h>
#include "declarations.h"

ArgumentParser::ParseResult parseCommandLine(HypervariableReconstructionArguments& arguments,
								    		 HypervariableReconstructionOptions& options,
								    		 int argc, char** argv) {
	
	ArgumentParser parser("ReconstructCDRs");

	/* Three required arguments */
	addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "[readFileNameL]"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "[readFileNameR]"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "[bestIsotypesFileName]"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "[organism]"));

	/* Options */
	addOption(parser, ArgParseOption("o", "output-file-name-prefix",
				 					 "prefix for output fasta file with reconstructed CDRs",
				 					  ArgParseArgument::STRING, "STRING"));
	addOption(parser, ArgParseOption("M", "maximum-reads-aligned",
									 "Maximum number of reads aligned before returning consensus (-1 means use all reads)",
									  ArgParseArgument::INTEGER, "INT"));
	addOption(parser, ArgParseOption("Hr", "max-h-mutation-rate",
									 "Maximum allowed mutation rate for hypervariable regions (CDR1 and CDR2)",
									 ArgParseArgument::DOUBLE, "DOUBLE"));

	addOption(parser, ArgParseOption("Fr", "max-f-mutation-rate",
									 "Maximum allowed mutation rate for IGV framework regions",
									 ArgParseArgument::DOUBLE, "DOUBLE"));

	addOption(parser, ArgParseOption("Hl", "min-h-overlap-length",
									 "Minimum overlap length for putative alignments between hypervariable regions (CDR1 and CDR2) and reads",
									 ArgParseArgument::INTEGER, "INT"));

	addOption(parser, ArgParseOption("Fl", "min-f-overlap-length",
									 "Minimum overlap length for putative alignments between IGV framework regions and reads",
									 ArgParseArgument::INTEGER, "INT"));

	// addOption(parser, ArgParseOption("V", "--verbose",
	// 								 "Show alignments in real time", ArgParseArgument))

	/* Defaults */
	setDefaultValue(parser, "output-file-name-prefix", "reconstructions");
	setDefaultValue(parser, "maximum-reads-aligned", -1);
	setDefaultValue(parser, "max-h-mutation-rate", 0.35);
	setDefaultValue(parser, "max-f-mutation-rate", 0.2);
	setDefaultValue(parser, "min-h-overlap-length", 5);
	setDefaultValue(parser, "min-f-overlap-length", 5);

	/* Parse command line */
	ArgumentParser::ParseResult res = parse(parser, argc, argv);

	if (res != ArgumentParser::PARSE_OK) {
        return res;
	}

	/* Get option values */
	getOptionValue(options.HMaxMutationRate, parser, "max-h-mutation-rate");
	getOptionValue(options.FMaxMutationRate, parser, "max-f-mutation-rate");
	getOptionValue(options.HMinOverlapLength, parser, "min-h-overlap-length");
	getOptionValue(options.FMinOverlapLength, parser, "min-f-overlap-length");
	getOptionValue(options.maximumReadsAligned, parser, "maximum-reads-aligned");
	getOptionValue(options.outputFileName, parser, "output-file-name-prefix");

	/* Get argument values */
	getArgumentValue(arguments.readFileNameL, parser, 0);
	getArgumentValue(arguments.readFileNameR, parser, 1);
	getArgumentValue(arguments.bestIsosFileName, parser, 2);
	getArgumentValue(arguments.organism, parser, 3);

	if (options.HMinOverlapLength <= 0 || options.FMinOverlapLength <= 0) {
		std::cerr << "ReconstructCDRs: You cannot use an overlap length option that is <= 0!\n";
		return ArgumentParser::PARSE_ERROR;
	} else if (options.HMaxMutationRate < 0 || options.FMaxMutationRate < 0) {
		std::cerr << "ReconstructCDRs: You cannot use a negative mutation rate! That doesn't make any sense...\n";
		return ArgumentParser::PARSE_ERROR;
	} else if (options.HMaxMutationRate >= 1 || options.FMaxMutationRate >= 1) {
		std::cerr << "ReconstructCDRs: Mutation rates greater than or eqaul to 1 are invalid\n";
		return ArgumentParser::PARSE_ERROR;
	} 

	if (!(arguments.organism == "human" || arguments.organism == "mouse")) {
		std::cerr << "ReconstructCDRs: " << arguments.organism << " is not a valid organism. Choose human or mouse\n";
		return ArgumentParser::PARSE_ERROR;
	}

	std::cout << std::flush;
	
	return ArgumentParser::PARSE_OK;
}

#endif