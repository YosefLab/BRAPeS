
#ifndef DECLARATIONS
#define DECLARATIONS

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/modifier.h>
#include <seqan/align.h>

#include <vector>
#include <map>
#include <regex>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

using namespace seqan;

static const unsigned int GAP_ORD_VAL = 5;
static const unsigned int MAX_REGEX_LENGTH = 1024;

typedef Dna5String TString;
typedef StringSet<TString> TStringSet;
typedef StringSet<CharString> TNameSet;
typedef Align<TString, ArrayGaps> TAlign;
typedef FragmentStore<> TStore;
typedef Row<TAlign>::Type TRow;
// typedef ModifiedString<ModifiedString<TString, ModComplementDna5>, ModReverse> ReverseComplement;

typedef std::map<CharString, TString> FragmentHashTable;
typedef std::map<CharString, int> RelativeReadCoordinateHashTable;

typedef struct HypervariableReconstructionOptions {
	double HMaxMutationRate;
	double FMaxMutationRate;
	unsigned int HMinOverlapLength;
	unsigned int FMinOverlapLength;
	int maximumReadsAligned;
	CharString outputFileName;
	bool verbose;
} HypervariableReconstructionOptions;

typedef struct HypervariableReconstructionArguments {
	CharString readFileNameL;
	CharString readFileNameR;
	CharString bestIsosFileName;
	CharString organism;
} HypervariableReconstructionArguments;

#endif