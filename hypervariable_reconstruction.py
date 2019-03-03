import os
import subprocess as sp  
from Bio import SeqIO
import pysam
import sys
import glob
import shutil
from Bio import pairwise2

CDR1_START = 78
CDR1_END = 114
CDR2_START = 165
CDR2_END = 195


reconstruction_path, best_iso, reads_file_1, reads_file_2, organism, reads_info, framework_extension, threshold_score, min_reads, max_reads, reads_before_moving_on, output_path_prefix, rsem, imgt_file = sys.argv[1:]

def removeFiles(path_plus_wildcards):
	for file in glob.glob(path_plus_wildcards):
		if os.path.isdir(file):
			shutil.rmtree(file)
		else:
			os.remove(file)

def genomicAADict(conserved_aa_file):
	return_dict = {}
	with open(conserved_aa_file, "r") as file:
		for line in file.readlines():
			key, value = line.rstrip().split()
			return_dict[key] = value
	return return_dict

def genomicNTDict(fasta_file):
	return_dict = {}
	for record in SeqIO.parse(fasta_file, "fasta"):
		key = record.id
		value = str(record.seq)
		return_dict[key] = value
	return return_dict

def IMGTdict(imgt_file):
	IMGT_dict = {}
	for record in SeqIO.parse(imgt_file, "fasta"):
		gene, version = record.id.split("|")[1].split("*")
		if version == "01":
			IMGT_dict[gene.replace("/", "")] = str(record.seq)
	return IMGT_dict

def getGenomicCDRs(imgt_sequence, framework_extension):
	CDR1 = imgt_sequence[CDR1_START - framework_extension:CDR1_END + framework_extension].replace(".", "").upper()
	CDR2 = imgt_sequence[CDR2_START - framework_extension:CDR2_END + framework_extension].replace(".", "").upper()
	return CDR1, CDR2

def makeSPCall(command_list, writeout=False):
	if not writeout:
		with open(os.devnull, 'w') as devnull:
			sp.call(command_list, stdout=devnull, stderr=devnull)
	else:
		sp.call(command_list)

def writeMessage(message):
	sys.stdout.write(message)
	sys.stdout.flush()

def runPipe(best_iso, reads_file_1, reads_file_2, imgt_file, framework_extension, threshold_score, min_reads, max_reads, reads_before_moving_on, output_path_prefix, reads_info, rsem):

	imgt_dict = IMGTdict(imgt_file)

	if not os.path.exists("ref") and not os.path.isdir("ref"):
		os.mkdir("ref")

	HVR_records = []
	reconstructed_iso_records = []

	for transcript_record in SeqIO.parse(best_iso, "fasta"):

		V_gene_symbol = transcript_record.id.split(".")[0]
		genomicCDR1 = ""
		genomicCDR2 = ""

		try:
			genomicCDR1, genomicCDR2 = getGenomicCDRs(imgt_dict[V_gene_symbol], framework_extension)
			#writeMessage("{} Genomic CDR1: {}\n".format(V_gene_symbol, genomicCDR1))
			#writeMessage("{} Genomic CDR2: {}\n".format(V_gene_symbol, genomicCDR2))
		except KeyError:
			writeMessage("\nSkipping... {} is not a {} gene\n".format(V_gene_symbol, organism))
			continue

		### Make individual transcript fasta file ###
		SeqIO.write(transcript_record, transcript_record.id + ".fa", "fasta")

		### Make rsem reference for transcript ###
		writeMessage("\nMaking RSEM reference for {}\n".format(transcript_record.id))
		makeSPCall([rsem + "rsem-prepare-reference", "--bowtie2", "-q", transcript_record.id + ".fa", "ref/" + transcript_record.id + ".ref"])

		### Calculate rsem expression results of paired-end reads on transcript ###
		writeMessage("Calculating expression for {}\n\n".format(transcript_record.id))
		makeSPCall([rsem + "rsem-calculate-expression", "-q", "-p", "4", "--no-qualities", "--bowtie2", "--paired-end",
				  reads_file_1, reads_file_2, "ref/" + transcript_record.id + ".ref", transcript_record.id])

		### Clean ###
		removeFiles("*results")
		removeFiles("*stat")


		reconstructedCDR1 = ""
		reconstructedCDR2 = ""
		
		writeMessage("\nAligning Reads...\n")
		makeSPCall([reconstruction_path + "build/SeqanReconstruct", transcript_record.id + ".transcript.bam", genomicCDR1, genomicCDR2, threshold_score, 
											min_reads, max_reads, reads_before_moving_on, str(transcript_record.seq), reads_info, str(framework_extension)])
		for record in SeqIO.parse("hvr.reconstructions.temp.fasta", "fasta"):
			record.id = transcript_record.id
			#writeMessage("{}    {}\n".format(record.id, record.seq))
			# if framework_extension != 0:
			# 	record.seq = record.seq[framework_extension:-framework_extension]
				
			if record.name == "CDR1":
				reconstructedCDR1 = record.seq
			elif record.name == "CDR2":
				reconstructedCDR2 = record.seq
			HVR_records.append(record)
		for record in SeqIO.parse("bcr.reconstructions.temp.fasta", "fasta"):
			record.id = transcript_record.id
			reconstructed_iso_records.append(record)

		removeFiles(transcript_record.id + "*")

	### clean ###
	removeFiles("*reconstructions.temp*")
	removeFiles("ref")

	SeqIO.write(HVR_records, output_path_prefix + ".hypervariable.regions.fasta", "fasta");
	SeqIO.write(reconstructed_iso_records, output_path_prefix + ".CDR1.CDR2.reconstructions.fasta", "fasta")

	# writeMessage("\nRanking BCR reconstructions...\n")
	# makeSPCall(["rsem-prepare-reference", "--bowtie2", "-q", output_path_prefix + ".BCR.reconstructions.fasta", "ref/" + output_path_prefix + ".ref"])
	# makeSPCall(["rsem-calculate-expression", "-q", "-p", "4", "--no-qualities", "--bowtie2", "--paired-end", 
	# 			reads_file_1, reads_file_2, "ref/" + output_path_prefix + ".ref", output_path_prefix])



if __name__ == "__main__":
	if not reconstruction_path.endswith("/"):
		reconstruction_path += "/"
	runPipe(best_iso, reads_file_1, reads_file_2, imgt_file, int(framework_extension), threshold_score, min_reads, max_reads, reads_before_moving_on, output_path_prefix, reads_info, rsem)











