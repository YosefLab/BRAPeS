#! /usr/bin/env python
import sys
import shutil
import os
import random
import argparse
import subprocess
import datetime
from Bio import SeqIO
from Bio import pairwise2
import pysam
import operator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def runBCRpipe(genome, output, bam, unmapped, bases, strand, numIterations,thresholdScore, minOverlap, rsem, bowtie2,
               path, sumF, NolowQ, samtools, top, byExp, readOverlap, oneSide, downsample,
               Hminus, Kminus, Lminus, skipHVR, Hr, Fr):
    checkParameters(strand, path, sumF, Hr, Fr)
    if path == './':
        path = os.getcwd()
    if not path.endswith('/'):
        path = path + '/'
    finalStatDict = dict()
    tcrFout = open(sumF + '.BCRs.txt','w')
    opened = False
    currFolder = os.path.abspath(os.path.dirname(sys.argv[0])) + '/'
    hvr_path = currFolder + 'HVR_recon'
    (fasta, bed, mapping, aaF) = getGenomeFiles(currFolder, genome)
    for cellFolder in os.listdir(path):
        fullPath = path + cellFolder + '/'
        if((os.path.exists(fullPath)) & (os.path.isdir(fullPath))):
            sys.stdout.write(str(datetime.datetime.now()) + " Working on: " + cellFolder + '\n')
            sys.stdout.flush()
            (foundBam, foundUnmapped, nbam, nunmapped, noutput) = formatFiles(fullPath, bam, unmapped, output)
            if not foundBam:
                sys.stderr.write(str(datetime.datetime.now()) + " There is not a bam in "
                                                                "this folder, moving to the next folder\n")
                sys.stderr.flush()
            elif not foundUnmapped:
                sys.stderr.write(str(datetime.datetime.now()) + " There is not an unmapped bam in "
                                                                "this folder, moving to the next folder\n")
                sys.stderr.flush()

            else:
                reconstruction = currFolder + '/vdj.alignment.bcr'
                runSingleCell(fasta, bed, noutput, nbam, nunmapped, mapping, bases, strand, reconstruction, aaF , numIterations, thresholdScore,
                            minOverlap, rsem, bowtie2, NolowQ, samtools, top, byExp, readOverlap, oneSide, downsample, genome,
                              Hminus, Kminus, Lminus, skipHVR, Hr, Fr,
                              hvr_path, genome)
                opened = addCellToTCRsum(cellFolder, noutput, opened, tcrFout)
                finalStatDict = addToStatDict(noutput, cellFolder, finalStatDict)
    sumFout = open(sumF + '.summary.txt','w')
    sumFout.write('sample\theavy\tkappa\tlambda\n')
    for cell in sorted(finalStatDict):
        fout = cell + '\t' + finalStatDict[cell]['heavy'] + '\t' + finalStatDict[cell]['kappa'] + '\t' + finalStatDict[cell]['lambda'] + '\n'
        sumFout.write(fout)
    sumFout.close()


def getGenomeFiles(currFolder, genome):
    dataFold = currFolder + 'Data/' + genome + '/'
    fasta = 'NA'
    bed = 'NA'
    mapping = 'NA'
    aaF = 'NA'
    if os.path.isdir(dataFold):
        for ff in os.listdir(dataFold):
            if not ff.startswith('.'):
                if ff.endswith('BCR.bed'):
                    bed = dataFold + ff
                elif ff.endswith('conserved.AA.txt'):
                    aaF = dataFold + ff
                elif ff.endswith('BCR.fa'):
                    fasta = dataFold + ff
                elif ff.endswith('gene.id.mapping.BCR.txt'):
                    mapping = dataFold + ff
        if (fasta == 'NA'):
            sys.exit(str(datetime.datetime.now()) + " Error! BCR fasta file is missing in the Data/genome folder\n")
        if (mapping == 'NA'):
            sys.exit(str(datetime.datetime.now()) + " Error! BCR gene id mapping file is missing in the Data/genome folder\n")
        if (bed == 'NA'):
            sys.exit(str(datetime.datetime.now()) + " Error! BCR bed file is missing in the Data/genome folder\n")
        if (aaF == 'NA'):
            sys.exit(str(datetime.datetime.now()) + " Error! BCR conserved AA file is missing in the Data/genome folder\n")
    else:
        sys.exit(str(datetime.datetime.now()) + " Error! Genome parameter is invalid, no folder named: " + dataFold + '\n')

    return(fasta, bed, mapping, aaF)

def addCellToTCRsum(cellFolder, noutput, opened, tcrFout):
    if os.path.isfile(noutput + '.summary.txt'):
        currOut = open(noutput + '.summary.txt','r')
        if not opened:
            opened = True
            head = currOut.readline()
            head = 'cell\t' + head
            tcrFout.write(head)
        else:
            currOut.readline()
        l = currOut.readline()
        while l != '':
            newL = cellFolder + '\t' + l
            tcrFout.write(newL)
            l = currOut.readline()
        currOut.close()
    return opened

def addToStatDict(noutput, cellFolder, finalStatDict):
    if cellFolder in finalStatDict:
        print("Error! %s appear more than once in final stat dictionary" % cellFolder)
    finalStatDict[cellFolder] = {'heavy':'Failed - found V and J segments but wasn\'t able to extend them',
                                 'kappa':'Failed - found V and J segments but wasn\'t able to extend them',
                                 'lambda':'Failed - found V and J segments but wasn\'t able to extend them'}
    if os.path.isfile(noutput + '.summary.txt'):
        currOut = open(noutput + '.summary.txt','r')
        msgK = 'None'
        msgH = 'None'
        msgL = 'None'
        currOut.readline()
        l = currOut.readline()
        while l != '':
            lArr = l.strip('\n').split('\t')
            chain = lArr[0]
            stat = lArr[1]
            if stat == 'Productive':
                if chain == 'heavy':
                    msgH = 'Productive'
                elif chain == 'kappa':
                    msgK = 'Productive'
                else:
                    msgL = 'Productive'
            elif stat == 'Productive (no 118 PHE/TRP found)':
                if chain == 'heavy':
                    msgH = 'Productive (no 118 PHE/TRP found)'
                elif chain == 'kappa':
                    msgK = 'Productive (no 118 PHE/TRP found)'
                else:
                    msgL = 'Productive (no 118 PHE/TRP found)'
            elif stat.startswith('Unproductive'):
                if chain == 'heavy':
                    if msgH != 'Productive':
                        msgH = 'Unproductive'
                elif chain == 'kappa':
                    if msgK != 'Productive':
                        msgK = 'Unproductive'
                else:
                    if msgL != 'Productive':
                        msgL = 'Unproductive'
            elif stat.startswith('Failed reconstruction'):
                if stat == 'Failed reconstruction - reached maximum number of iterations':
                    if chain == 'heavy':
                        if msgH == 'None':
                            msgH = 'Failed - reconstruction didn\'t converge'
                    elif chain == 'kappa':
                        if msgK == 'None':
                            msgK = 'Failed - reconstruction didn\'t converge'
                    else:
                        if msgL == 'None':
                            msgL = 'Failed - reconstruction didn\'t converge'
                elif stat == 'Failed reconstruction - V and J segment do not overlap':
                    if chain == 'heavy':
                        if msgH == 'None':
                            msgH = 'Failed - V and J reconstruction don\'t overlap'
                    elif chain == 'kappa':
                        if msgK == 'None':
                            msgK = 'Failed - V and J reconstruction don\'t overlap'
                    else:
                        if msgL == 'None':
                            msgL = 'Failed - V and J reconstruction don\'t overlap'
                            msgL = 'Failed - V and J reconstruction don\'t overlap'

            l = currOut.readline()
        currOut.close()
        if msgH == 'None':
            heavyJunc = noutput + '.heavy.junctions.txt'
            if (os.path.isfile(heavyJunc) == True):
                if os.stat(heavyJunc).st_size == 0:
                    msgH = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msgH = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msgH = 'Failed - didn\'t find any V and J segments in original mapping'
        if msgK == 'None':
            kappaJunc = noutput + '.kappa.junctions.txt'
            if (os.path.isfile(kappaJunc) == True):
                if os.stat(kappaJunc).st_size == 0:
                    msgK = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msgK = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msgK = 'Failed - didn\'t find any V and J segments in original mapping'
        if msgL == 'None':
            lambdaJunc = noutput + '.lambda.junctions.txt'
            if (os.path.isfile(lambdaJunc) == True):
                if os.stat(lambdaJunc).st_size == 0:
                    msgL = 'Failed - didn\'t find any V and J segments in original mapping'
                else:
                    msgL = 'Failed - found V and J segments but wasn\'t able to extend them'
            else:
                msgL = 'Failed - didn\'t find any V and J segments in original mapping'

    else:
        kappaJunc = noutput + '.kappa.junctions.txt'
        heavyJunc = noutput + '.heavy.junctions.txt'
        lambdaJunc = noutput + '.lambda.junctions.txt'
        if os.path.isfile(kappaJunc) == True:
            if os.stat(kappaJunc).st_size == 0:
                msgK = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgK = 'Failed - didn\'t find any V and J segments in original mapping'
        if (os.path.isfile(heavyJunc) == True):
            if os.stat(heavyJunc).st_size == 0:
                msgH = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgH = 'Failed - didn\'t find any V and J segments in original mapping'
        if (os.path.isfile(lambdaJunc) == True):
            if os.stat(lambdaJunc).st_size == 0:
                msgL = 'Failed - didn\'t find any V and J segments in original mapping'
        else:
            msgL = 'Failed - didn\'t find any V and J segments in original mapping'

    finalStatDict[cellFolder]['heavy'] = msgH
    finalStatDict[cellFolder]['kappa'] = msgK
    finalStatDict[cellFolder]['lambda'] = msgL
    return finalStatDict


def formatFiles(fullPath, bam, unmapped, output):
    foundBam = True
    foundUnmapped = True
    nbam = fullPath + bam
    if bam.startswith('/'):
        nbam = fullPath + bam[1:]
    elif bam.startswith('./'):
        nbam = fullPath + bam[2:]
    nunmapped = fullPath + unmapped
    if unmapped.startswith('/'):
        nunmapped = fullPath + unmapped[1:]
    if unmapped.startswith('./'):
        nunmapped = fullPath + unmapped[2:]
    if ((os.path.isfile(nunmapped)) & (os.path.isfile(nbam))):
        noutput = makeOutputDir(output, fullPath)
    else:
        if (os.path.isfile(unmapped)):
            foundUnmapped = False
        else:
            foundBam = False
        noutput = output
    return (foundBam, foundUnmapped, nbam, nunmapped, noutput)

def makeOutputDir(output, fullPath):
    noutput = output
    if output.startswith('/'):
        noutput = output[1:]
    if output.startswith('./'):
        noutput = output[2:]
    if output.endswith('/'):
        noutput = output[:-1]
    if output.find('/') != -1:
        outArr = noutput.split('/')
        currPath = fullPath
        for i in range(0,len(outArr)-1):
            currPath = currPath + outArr[i] + '/'
            if not os.path.exists(currPath):
                os.makedirs(currPath)
    noutput = fullPath + noutput
    return noutput


def runSingleCell(fasta, bed, output, bam, unmapped, mapping, bases, strand, reconstruction, aaF , numIterations, thresholdScore, minOverlap,
                  rsem, bowtie2, NolowQ, samtools, top, byExp, readOverlap, oneSide, downsample, organism,
                  Hminus, Kminus, Lminus, skipHVR, Hr, Fr,
                  hvr_path, genome):
    idNameDict = makeIdNameDict(mapping)
    fastaDict = makeFastaDict(fasta)
    vdjDict = makeVDJBedDict(bed, idNameDict)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing heavy chain\n")
    sys.stdout.flush()
    unDictHeavy = analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 'H', strand, NolowQ, top, byExp, readOverlap, downsample, organism, Hminus, Kminus, Lminus)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing kappa chain\n")
    sys.stdout.flush()
    unDictKappa = analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 'K', strand, NolowQ, top, byExp, readOverlap, downsample, organism, Hminus, Kminus, Lminus)
    sys.stdout.write(str(datetime.datetime.now()) + " Pre-processing lambda chain\n")
    sys.stdout.flush()
    unDictLambda = analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, 'L', strand, NolowQ, top, byExp, readOverlap, downsample, organism, Hminus, Kminus, Lminus)
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing heavy chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.heavy.mapped.and.unmapped.fa', output + '.heavy.junctions.txt', output + '.reconstructed.junctions.heavy.fa', str(numIterations), str(thresholdScore), str(minOverlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing kappa chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.kappa.mapped.and.unmapped.fa', output + '.kappa.junctions.txt',
        output + '.reconstructed.junctions.kappa.fa', str(numIterations), str(thresholdScore), str(minOverlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Reconstructing lambda chains\n")
    sys.stdout.flush()
    subprocess.call([reconstruction, output + '.lambda.mapped.and.unmapped.fa', output + '.lambda.junctions.txt',
        output + '.reconstructed.junctions.lambda.fa', str(numIterations), str(thresholdScore), str(minOverlap)])
    sys.stdout.write(str(datetime.datetime.now()) + " Creating full BCR sequences\n")
    sys.stdout.flush()
    fullTcrFileHeavy = output + '.heavy.full.BCRs.fa'
    tcrF = output + '.reconstructed.junctions.heavy.fa'
    (cSeq, cName, cId) = findBestCtoAppend(vdjDict['heavy']['C'], idNameDict, fastaDict, bam, downsample)
    createTCRFullOutput(fastaDict, tcrF, fullTcrFileHeavy, bases, idNameDict, cSeq, cName, cId, oneSide)
    fullTcrFileKappa = output + '.kappa.full.BCRs.fa'
    tcrF = output + '.reconstructed.junctions.kappa.fa'
    (cSeq, cName, cId) = findBestCtoAppend(vdjDict['kappa']['C'], idNameDict, fastaDict, bam, downsample)
    createTCRFullOutput(fastaDict, tcrF, fullTcrFileKappa , bases, idNameDict, cSeq, cName, cId, oneSide)
    fullTcrFileLambda = output + '.lambda.full.BCRs.fa'
    tcrF = output + '.reconstructed.junctions.lambda.fa'
    (cSeq, cName, cId) = findBestCtoAppend(vdjDict['lambda']['C'], idNameDict, fastaDict, bam, downsample)
    createTCRFullOutput(fastaDict, tcrF, fullTcrFileLambda , bases, idNameDict, cSeq, cName, cId, oneSide)
    sys.stdout.write(str(datetime.datetime.now()) + " Running RSEM to quantify expression of all possible isotypes\n")
    sys.stdout.flush()
    outDirInd = output.rfind('/')
    if outDirInd != -1:
        outDir = output[:outDirInd+1]
    else:
        outDir = os.getcwd()
    runRsem(outDir, rsem, bowtie2, fullTcrFileHeavy, fullTcrFileKappa,fullTcrFileLambda, output, samtools, False)
    pickFinalIsoforms(fullTcrFileHeavy, fullTcrFileKappa, fullTcrFileLambda , output)
    bestHeavy = output + '.heavy.full.BCRs.bestIso.fa'
    bestKappa = output + '.kappa.full.BCRs.bestIso.fa'
    bestLambda = output + '.lambda.full.BCRs.bestIso.fa'
    sys.stdout.write(str(datetime.datetime.now()) + " Finding productive CDR3\n")
    sys.stdout.flush()
    aaDict = makeAADict(aaF)
    if os.path.isfile(bestHeavy):
        fDictHeavy = findCDR3(bestHeavy, aaDict, fastaDict, False )
    else:
        fDictHeavy = dict()
    if os.path.isfile(bestKappa):
        fDictKappa = findCDR3(bestKappa, aaDict, fastaDict, False )
    else:
        fDictKappa = dict()
    if os.path.isfile(bestLambda):
        fDictLambda = findCDR3(bestLambda, aaDict, fastaDict, False )
    else:
        fDictLambda = dict()
    if skipHVR:
        sys.stdout.write(str(datetime.datetime.now()) + " Skipping CDR1 and CDR2 reconstruction\n")
        sys.stdout.flush()
    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Performing CDR1 and CDR2 reconstruction\n")
        sys.stdout.flush()
        if os.path.isfile(bestHeavy):
            bestHeavy = runHVRrecon(bestHeavy, 'heavy', hvr_path, output, genome, downsample, Hr, Fr)
        if os.path.isfile(bestKappa):
            bestKappa = runHVRrecon(bestKappa, 'kappa',hvr_path, output, genome, downsample, Hr, Fr)
        if os.path.isfile(bestLambda):
            bestLambda = runHVRrecon(bestLambda, 'lambda', hvr_path, output, genome, downsample, Hr, Fr)
        runRsem(outDir, rsem, bowtie2, bestHeavy, bestKappa,bestLambda, output, samtools, True)
    if skipHVR:
        heavyRsemOut = output + '.heavy.rsem.out.genes.results'
        kappaRsemOut = output + '.kappa.rsem.out.genes.results'
        lambdaRsemOut = output + '.lambda.rsem.out.genes.results'
        heavyBam = output + '.heavy.rsem.out.transcript.sorted.bam'
        kappaBam = output + '.kappa.rsem.out.transcript.sorted.bam'
        lambdaBam = output + '.lambda.rsem.out.transcript.sorted.bam'
    else:
        if os.path.isfile(bestHeavy):
            fDictHeavy = findCDR3(bestHeavy, aaDict, fastaDict, True )
        else:
            fDictHeavy = dict()
        if os.path.isfile(bestKappa):
            fDictKappa = findCDR3(bestKappa, aaDict, fastaDict, True )
        else:
            fDictKappa = dict()
        if os.path.isfile(bestLambda):
            fDictLambda = findCDR3(bestLambda, aaDict, fastaDict, True )
        else:
            fDictLambda = dict()
        heavyRsemOut = output + '.afterCDRrec.heavy.rsem.out.genes.results'
        kappaRsemOut = output + '.afterCDRrec.kappa.rsem.out.genes.results'
        lambdaRsemOut = output + '.afterCDRrec.lambda.rsem.out.genes.results'
        heavyBam = output + '.afterCDRrec.heavy.rsem.out.transcript.sorted.bam'
        kappaBam = output + '.afterCDRrec.kappa.rsem.out.transcript.sorted.bam'
        lambdaBam = output + '.afterCDRrec.lambda.rsem.out.transcript.sorted.bam'
    sys.stdout.write(str(datetime.datetime.now()) + " Writing results to summary file\n")
    sys.stdout.flush()
    makeSingleCellOutputFile(fDictHeavy, fDictKappa, fDictLambda, output, heavyRsemOut,kappaRsemOut , lambdaRsemOut,
                             heavyBam, kappaBam, lambdaBam, fastaDict, unDictHeavy,unDictKappa, unDictLambda, idNameDict,
                             downsample)



def runHVRrecon(bestIso, chain, hvr_path, output, genome, downsample, Hr, Fr):
    pref = bestIso.replace('.fa','')
    if genome.startswith('hg'):
        org = 'human'
    elif genome.startswith('mm'):
        org = 'mouse'
    else:
        print("genome is not human or mouse, not running any CDR1/2 reconstruction")
    if not downsample:
        subprocess.call([hvr_path + '/build/ReconstructCDRs','-Hr',str(Hr),'-Fr',str(Fr),'-o',
                         pref + '.CDR1.CDR2.reconstructions.fasta',output + '.' + chain + '.R1.fa',
                         output + '.' + chain + '.R2.fa',bestIso, org])
    else:
        subprocess.call([hvr_path + '/build/ReconstructCDRs','-Hr',str(Hr),'-Fr',str(Fr),'-M','40000','-o',
                         pref + '.CDR1.CDR2.reconstructions.fasta',output + '.' + chain + '.R1.fa',
                         output + '.' + chain + '.R2.fa',bestIso, org])
    curBestIso = pref + '.CDR1.CDR2.reconstructions.fasta.BCRs.fasta'
    return curBestIso




def findBestCtoAppend(cArr, idNameDict, fastaDict, mappedBam, downsample):
    bestC = 0
    countNum = 0
    for i in range(0,len(cArr)):
        lArr = cArr[i].strip('\n').split('\t')
        chr = lArr[0]
        start = int(lArr[1])
        end = int(lArr[2])
        curCounts = findCountsInRegion(mappedBam, start, end, chr, downsample)
        if curCounts > countNum:
            countNum = curCounts
            bestC = i
    return getCInfo(cArr[bestC],idNameDict, fastaDict)



def makeSingleCellOutputFile(heavyDict, kappaDict, lambdaDict, output, heavyRsem, kappaRsem , lambdaRsem,
                             heavyBam, kappaBam, lambdaBam, fastaDict, unDictHeavy, unDictKappa, unDictLambda, idNameDict,
                             downsample):
    outF = open(output + '.summary.txt', 'w')
    outF.write('Chain\tStatus\tRank of BCR\tV\tJ\tC\tCDR3 NT\tCDR3 AA\t#reads in BCR\t#reads in CDR3\t#reads in V\t#reads in J\t#reads in C\t%unmapped reads used in the reconstruction\t# unmapped reads used in the reconstruction\t%unmapped reads in CDR3\t#unmapped reads in CDR3\tV ID\tJ ID\tC ID\n')
    if (len(heavyDict) > 0):
        writeChain(outF, 'heavy',heavyDict,heavyRsem, heavyBam, fastaDict,unDictHeavy, output, idNameDict, downsample)
    if (len(kappaDict) > 0):
        writeChain(outF,'kappa',kappaDict, kappaRsem, kappaBam, fastaDict, unDictKappa, output, idNameDict, downsample)
    if (len(lambdaDict) > 0):
        writeChain(outF,'lambda',lambdaDict, lambdaRsem, lambdaBam, fastaDict, unDictLambda, output, idNameDict, downsample)
    outF.close()


def writeChain(outF, chain,cdrDict,rsemF, bamF, fastaDict, unDict, output, idNameDict, downsample):
    writtenArr = []
    if os.path.exists(rsemF):
        noRsem = False
        (rsemDict,unRsemDict) = makeRsemDict(rsemF, cdrDict)
    else:
        noRsem = True
    for tcr in cdrDict:
        jStart = -1
        cdrInd = -1
        cInd = -1
        if cdrDict[tcr]['stat'].startswith('Productive'):
            isProd = True
        else:
            isProd = False
        fLine = chain + '\t' + cdrDict[tcr]['stat'] + '\t'
        if noRsem:
            rank = 'NA'
        else:
            #print rsemDict
            #print unRsemDict
            rank = getRank(tcr, rsemDict, unRsemDict, isProd, noRsem)
        fLine += str(rank) + '\t'
        nameArr = tcr.split('.')
        fLine += nameArr[0] + '\t' + nameArr[1] + '\t' + nameArr[2] + '\t'
        fLine += cdrDict[tcr]['CDR3 NT'] + '\t' + cdrDict[tcr]['CDR3 AA'] + '\t'
        fullSeq = cdrDict[tcr]['Full Seq'].upper()
        if not noRsem:
            totalCount = findCountsInRegion(bamF, 0, len(fullSeq), tcr, downsample)
            fLine += str(totalCount) + '\t'
            cName = nameArr[5]
            while cName.endswith('_2'):
                cName = cName[:-2]
            cSeq = fastaDict[cName].upper()
            cInd = fullSeq.rfind(cSeq)
            if cInd == -1:
                sys.stderr.write(str(datetime.datetime.now()) + ' Warning! could not find C segment sequence in the full sequence\n')
                sys.stderr.flush()
                cCounts = 'NA'
            else:
                cCounts = findCountsInRegion(bamF, cInd, len(fullSeq), tcr, downsample)
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                cdrInd = fullSeq.find(cdrDict[tcr]['CDR3 NT'].upper())
            else:
                cdrInd = -1
            if ((cdrInd == -1) & (cdrDict[tcr]['CDR3 NT'] != 'NA')):
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! Cound not find CDR3 NT sequence in the full sequence\n')
                sys.stderr.flush()
            if cdrInd != -1:
                cdrCounts = findCountsInRegion(bamF, cdrInd, cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr, downsample)
                jStart = cdrInd + len(cdrDict[tcr]['CDR3 NT'])
                if cInd != -1:
                    if cInd < jStart:
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! index for C segement start is before end of J segment index, ignoring C index position\n')
                        sys.stderr.flush()
                        jCounts = findCountsInRegion(bamF, jStart, jStart+50, tcr, downsample)
                    else:
                        jCounts = findCountsInRegion(bamF, jStart, cInd, tcr, downsample)
                else:
                    jCounts = findCountsInRegion(bamF, jStart, jStart+50, tcr, downsample)
                vCounts = findCountsInRegion(bamF, 0, cdrInd, tcr, downsample)
                fLine += str(cdrCounts) + '\t' + str(vCounts) + '\t' + str(jCounts) + '\t' + str(cCounts) + '\t'
            else:
                fLine += 'NA\tNA\tNA\t' + str(cCounts) + '\t'
            vId = nameArr[3]
            jId = nameArr[4]
            cId = nameArr[5]
            if cdrDict[tcr]['CDR3 NT'] != 'NA':
                (unDictRatioCDR, unCDRcount) = getUnDictRatio(bamF, cdrInd , cdrInd + len(cdrDict[tcr]['CDR3 NT']), tcr, unDict, downsample)
                (unDictRatioALL, unAllcount) = getUnDictRatio(bamF, 0 , len(fullSeq), tcr, unDict, downsample)
                fLine += str(unDictRatioALL) + '\t' + str(unAllcount) + '\t' + str(unDictRatioCDR) + '\t' + str(unCDRcount) + '\t'
            else:
                fLine += 'NA\tNA\tNA\tNA\t'
            writtenArr.append(vId)
            writtenArr.append(jId)
            fLine += vId + '\t' + jId + '\t' + cId + '\n'
        else:
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t' + nameArr[3] + '\t' + nameArr[4] + '\t' + nameArr[5] + '\n'
            #print fLine

        outF.write(str(fLine))
    writeFailedReconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict )

def writeFailedReconstructions(outF, chain, writtenArr, output, idNameDict, fastaDict):
    recF = output + '.reconstructed.junctions.' + chain + '.fa'
    if os.path.isfile(recF):
        f = open(recF, 'rU')
        segDict = dict()
        for tcrRecord in SeqIO.parse(f, 'fasta'):
            tcrSeq = str(tcrRecord.seq)
            if tcrSeq.find('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN') != -1:
                status = 'Failed reconstruction - reached maximum number of iterations'
                segDict = addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
            elif tcrSeq.find('NNNN') != -1:
                status = 'Failed reconstruction - V and J segment do not overlap'
                segDict = addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict)
        f.close()
        if len(segDict) > 0:
            writeSegDict(segDict, outF, chain)


def writeSegDict(segDict, outF, chain):
    for seg in segDict:
        currDict = segDict[seg]
        pairs = ''
        for pair in currDict['pairs']:
            pairs += pair + '.'
        pairs = pairs[:-1]
        if currDict['len'] > 0:
            fLine = chain + '\t' + currDict['status'] + '\t'
            rank = findCurrRank(segDict, seg, currDict['len'])
            fLine += str(rank) + '\t'
            if currDict['type'] == 'V':
                fLine += currDict['name'] + '\t' + 'paired with: ' + pairs + '\t'
            else:
                fLine +=  'paired with: ' + pairs + '\t' + currDict['name'] + '\t'
            fLine += 'NA\t' + currDict['seq'] + '\tNA\tNA\tNA\t'
            fLine += 'NA\tNA\tNA\tNA\tNA\tNA\tNA\t'
            if currDict['type'] == 'V':
                fLine += seg + '\tNA\tNA\n'
            else:
                fLine += 'NA\t' + seg + '\tNA\n'
            outF.write(str(fLine))


def findCurrRank(segDict, seg, currLen):
    rank = 1
    for s in segDict:
        if s != seg:
            if segDict[s]['len'] > currLen:
                rank += 1
    return rank

def addSegmentsToDict(segDict, status, writtenArr, tcrRecord, idNameDict, fastaDict):
    head = tcrRecord.id
    headArr = head.split('.')
    vId = headArr[0]
    jId = headArr[1].split('(')[0]
    currSeqArr = tcrRecord.seq.split('N')
    vSeq = currSeqArr[0]
    jSeq = currSeqArr[-1]
    minLen = min(len(fastaDict[vId]),len(fastaDict[jId]))
    tupArr = [(vId, vSeq),(jId, jSeq)]
    for i in range(0,len(tupArr)):
        (id,seq) = tupArr[i]
        if id not in writtenArr:
            if id in segDict:
                if i == 0:
                    if idNameDict[jId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[jId])
                    if str(segDict[id]['seq'][-20:]) != str(seq[-20:]):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same V-segment %s\n' % id)
                        sys.stderr.flush()
                else:
                    if idNameDict[vId] not in segDict[id]['pairs']:
                        segDict[id]['pairs'].append(idNameDict[vId])
                    if str(segDict[id]['seq'][:20]) != str(seq[:20]):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconstructed two different sequences from the same J-segment %s\n' % id)
                        sys.stderr.flush()
            else:
                segDict[id] = dict()
                segDict[id]['status'] = status
                segDict[id]['seq'] = seq
                segDict[id]['len'] = len(seq) - minLen
                segDict[id]['pairs'] = []

                if i == 0:
                    segDict[id]['type'] = 'V'
                    segDict[id]['pairs'].append(idNameDict[jId])
                else:
                    segDict[id]['type'] = 'J'
                    segDict[id]['pairs'].append(idNameDict[vId])
                segDict[id]['name'] = idNameDict[id]

    return segDict






def getRank(tcr, rsemDict, unRsemDict, isProd, noRsem):
    if isProd:
        currDict = rsemDict
    else:
        #currDict = unRsemDict
        return 'NA'
    if not noRsem:
        currCount = currDict[tcr]
        rank = 1
        for rec in currDict:
            if rec != tcr:
                if unRsemDict[rec] > currCount:
                    rank += 1
        return rank
    else:
        return 'NA'


def getUnDictRatio(bamF, start, end, tcr, unDict, downsample):
    unMappedCount = 0
    if downsample:
        return(0,unMappedCount)
    usedArr = []
    mappedFile = pysam.AlignmentFile(bamF,"rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1 :
            newName =  read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in usedArr:
            usedArr.append(newName)
            if newName in unDict:
                unMappedCount += 1
            #if downsample:
            #   if unMappedCount > 100:
            #        if len(unDict) == 0:
            #            return(0,unMappedCount)
            #       else:
            #            return (float(float(unMappedCount)/len(unDict)), unMappedCount)
    mappedFile.close()
    if len(unDict) == 0:
        return(0,unMappedCount)
    else:
        return (float(float(unMappedCount)/len(unDict)), unMappedCount)


def findCountsInRegion(bamF, start, end, tcr, downsample):
    readsArr = dict()
    mappedFile = pysam.AlignmentFile(bamF,"rb")
    readsIter = mappedFile.fetch(tcr, start, end)
    for read in readsIter:
        if read.is_read1 :
            newName =  read.query_name + '_1'
        else:
            newName = read.query_name + '_2'
        if newName not in readsArr:
            readsArr[newName] = '1'
        if downsample:
            if len(readsArr) > 100:
                mappedFile.close()
                return len(readsArr)
    mappedFile.close()
    counts = len(readsArr)
    return counts

def makeRsemDict(rsemF, cdrDict):
    fDict = dict()
    unDict = dict()
    f = open(rsemF,'r')
    f.readline()
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        name = lArr[1]
        if name in cdrDict:
            if cdrDict[name]['stat'].startswith('Productive'):
                fDict[name] = float(lArr[4])
            unDict[name] = float(lArr[4])
        l = f.readline()
    f.close()
    return (fDict,unDict)



def findCDR3(fasta, aaDict, vdjFaDict, roundtwo):
    f = open(fasta, 'rU')
    fDict = dict()
    for record in SeqIO.parse(f, 'fasta'):
        if record.id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! same name for two fasta entries %s\n' % record.id)
            sys.stderr.flush()
        else:
            idArr = record.id.split('.')
            vSeg = idArr[0]
            jSeg = idArr[1]
            if ((vSeg in aaDict) & (jSeg in aaDict)):
                currDict = findVandJaaMap(aaDict[vSeg],aaDict[jSeg],record.seq, roundtwo)
            else:
                if vSeg in aaDict:
                    newVseg = aaDict[vSeg]
                else:
                    vId = idArr[3]
                    currSeq = vdjFaDict[vId]
                    newVseg = getBestVaa(Seq(currSeq))
                if jSeg in aaDict:
                    newJseg = aaDict[jSeg]
                else:
                    jId = idArr[4]
                    currSeq= vdjFaDict[jId]
                    newJseg = getBestJaa(Seq(currSeq))
                currDict = findVandJaaMap(newVseg,newJseg,record.seq, roundtwo)
            fDict[record.id] = currDict
    f.close()
    return fDict


## Search for one of the following motifs in the first 8 AA of the J seg: FXG, GXG.
## If both dont exsits, simply report the position of an F AA, if one exsits.
def getBestJaa(currSeq):
    firstSeq = getNTseq(currSeq)
    secondSeq = getNTseq(currSeq[1:])
    thirdSeq = getNTseq(currSeq[2:])
    pos = 10
    seq = ''
    found = False
    for s in [firstSeq, secondSeq, thirdSeq]:
        tempSeq = s[:8]
        indF = tempSeq.find('F')
        if indF != -1:
            if ((indF+3) <= len(s)):
                if s[indF+3] == 'G':
                    found = True
                    if indF < pos:
                        pos = indF
                        seq = s[indF:]
        indG = tempSeq.find('G')
        if (indG != -1):
            if found == False:
                if ((indG + 2) <= len(s)):
                    if s[indG+2] == 'G':
                        if indG < pos:
                            found = True
                            seq = s[indG:]
                            pos = indG
        if ((found == False) & (indF != -1)):
            seq = s[indF:]
    if seq != '':
        return seq
    else:
        return firstSeq

def getBestVaa(currSeq):
    firstSeq = getNTseq(currSeq)
    secondSeq = getNTseq(currSeq[1:])
    thirdSeq = getNTseq(currSeq[2:])
    pos = 10
    seq = ''
    for s in [firstSeq, secondSeq, thirdSeq]:
        tempSeq = s[-8:]
        ind = tempSeq.find('C')
        stopInd = tempSeq.find('*')
        if ((ind != -1) & (stopInd == -1)):
            if ind < pos:
                goodRF = isGoodRF(s)
                if goodRF:
                    pos = ind
                    seq = s[:-8+ind + 1]
    if seq != '':
        return seq
    else:
        return firstSeq

def isGoodRF(s):
    mInd = s.find('M')
    if mInd == -1:
        return False
    stopInd = s.find('*')
    if stopInd == -1:
        return True
    stopIndNext = s[stopInd+1:].find('*')
    while stopIndNext != -1:
        stopInd = stopIndNext + stopInd + 1
        stopIndNext = s[stopInd+1:].find('*')
        mInd = s[stopInd+1:].find('M')
        mInd = mInd + stopInd + 1
    if mInd != -1:
        return True
    else:
        return False



def findVandJaaMap(vSeg,jSeg,fullSeq, roundtwo):
    fDict = dict()
    firstSeq = getNTseq(fullSeq)
    secondSeq = getNTseq(fullSeq[1:])
    thirdSeq = getNTseq(fullSeq[2:])
    ntArr = [fullSeq, fullSeq[1:],fullSeq[2:]]
    aaSeqsArr = [firstSeq, secondSeq, thirdSeq]
    cdrArr = []
    posArr = []
    fPharr = []
    for aaSeq in aaSeqsArr:
        #print "new aaseq"
        if not roundtwo:
            (cdr, pos, curPh, foundC) = getCDR3(aaSeq, vSeg,jSeg)
        else:
            (cdr, pos, curPh, foundC) = getCDR3secondPass(aaSeq, vSeg,jSeg)
        cdrArr.append(cdr)
        posArr.append(pos)
        fPharr.append(curPh)
    maxLen = 0
    bestCDR = ''
    bestSeq = ''
    hasStop = False
    bestPos = -1
    bestCDRnt = ''
    foundGood = False
    vPos = -1
    jPos = -1
    fPh = False
    for i in range(0,3):
        if posArr[i] != -1:
            if ((cdrArr[i] != 'Only J') & (cdrArr[i] != 'Only V')):
                if len(cdrArr[i]) > maxLen:
                    if cdrArr[i].find('*') == -1:
                        foundGood = True
                        bestCDR = cdrArr[i]
                        bestPos = posArr[i]
                        maxLen = len(cdrArr[i])
                        bestSeq = ntArr[i]
                        fPh = fPharr[i]
                    else:
                        if maxLen == 0:
                            foundGood = True
                            bestPos = posArr[i]
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
                            hasStop = True
                            fPh = fPharr[i]
                else:
                    if hasStop == True:
                        if cdrArr[i].find('*') == -1:
                            foundGood = True
                            bestPos = posArr[i]
                            hasStop = False
                            bestCDR = cdrArr[i]
                            maxLen = len(cdrArr[i])
                            bestSeq = ntArr[i]
                            fPh = fPharr[i]
            else:
                if not foundGood:
                    fPh = fPharr[i]
                    if (cdrArr[i] == 'Only J'):
                        jPos = posArr[i]-i
                    elif (cdrArr[i] == 'Only V'):
                        vPos = posArr[i]-i
    if ((vPos != -1) & (jPos != -1) & (not foundGood)):
        bestCDRnt = fullSeq[3*vPos:3*jPos]
        bestCDR = 'NA'
    elif bestPos != -1:
        bestCDRnt = bestSeq[3*bestPos : 3*bestPos+3*len(bestCDR)]
    if bestCDR.find('*') != -1:
        stat = 'Unproductive - stop codon'
    elif fPh:
        stat = 'Productive (no 118 PHE/TRP found)'
    else:
        stat = 'Productive'
    if maxLen == 0:
        if (('Only J' in cdrArr) & ('Only V' in cdrArr)):
            stat = 'Unproductive - Frame shift'
        else:
            if (('Only J' not in cdrArr) & ('Only V' in cdrArr)):
                stat = 'Unproductive - found only V segment'
            elif (('Only J' in cdrArr) & ('Only V' not in cdrArr)):
                if foundC == "notFoundButV":
                    stat = 'Productive (2nd-CYS 104 not identified)'
                else:
                    stat = 'Unproductive - found only J segment'
            elif (('Only J' not in cdrArr) & ('Only V' not in cdrArr)):
                stat = 'Unproductive - didn\'t find V and J segment'
            else:
                stat = 'Unproductive'
            bestCDR = 'NA'
            bestCDRnt = 'NA'
    fDict['stat'] = stat
    fDict['CDR3 AA'] = bestCDR
    fDict['CDR3 NT'] = bestCDRnt
    fDict['Full Seq'] = fullSeq
    return fDict


def getCDR3secondPass(aaSeq, vSeq, jSeq):
    pos = -1
    bestScore = 0
    foundC = "notFound"
    alignments = pairwise2.align.localms(aaSeq,vSeq,1,-2,-1,-.1)
    for align in alignments:
        score = align[2]
        if float(score) >= bestScore:
            bestScore = float(score)
            noGapPos = -1
            for i in range(0,len(align[1])):
                if align[1][i] != '-':
                    noGapPos = i
            counter = 0
            for i in range(0,noGapPos+1):
                if align[0][i] != '-':
                    counter += 1
            if pos < counter:
                pos = counter
            if align[0][noGapPos] == 'C':
                foundC = "Found"
            else:
                foundC = "notFoundButV"
    jPos = -1
    minDistJ = 5
    curPh = False
    for j in range(pos+1, len(aaSeq) - len(jSeq) + 1):
        subAA = aaSeq[j: j + len(jSeq)]
        if len(subAA) != len(jSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong subj length\n')
            sys.stderr.flush()
        dist = 0
        for m in range(0,len(jSeq)):
            if jSeq[m] != subAA[m]:
                dist += 1
        if (dist <= minDistJ):
            if isLegal(subAA):
                jPos = j
                minDistJ = dist
                curPh = False
            else:
                if dist < minDistJ:
                    curPh = True
                    jPos = j
                    minDistJ = dist
    if pos == -1:
        if jPos != -1:
            return('Only J', jPos, curPh, foundC)
        else:
            return('No V/J found', -1, curPh, foundC)
    else:
        if jPos == -1:
            return('Only V',pos, curPh, foundC)
    return(aaSeq[pos:jPos], pos, curPh, foundC)



def getCDR3(aaSeq, vSeq, jSeq):
    minDist = 24
    pos = -1
    bestScore = 0
    foundC = "notFound"
    alignments = pairwise2.align.localms(aaSeq,vSeq,1,-2,-1,-.1)
    for align in alignments:
        score = align[2]
        if float(score) > len(vSeq)-minDist:
            noGapPos = -1
            for i in range(0,len(align[1])):
                if align[1][i] != '-':
                    noGapPos = i
            counter = 0
            for i in range(0,noGapPos+1):
                if align[0][i] != '-':
                    counter += 1
            pos = counter
            if align[0][noGapPos] == 'C':
                foundC = "Found"
            else:
                foundC = "notFoundButV"
    jPos = -1
    minDistJ = 4
    curPh = False
    for j in range(pos+1, len(aaSeq) - len(jSeq) + 1):
        subAA = aaSeq[j: j + len(jSeq)]
        if len(subAA) != len(jSeq):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Wrong subj length\n')
            sys.stderr.flush()
        dist = 0
        for m in range(0,len(jSeq)):
            if jSeq[m] != subAA[m]:
                dist += 1
        if (dist <= minDistJ):
            if isLegal(subAA):
                jPos = j
                minDistJ = dist
                curPh = False
            else:
                if dist < minDistJ:
                    curPh = True
                    jPos = j
                    minDistJ = dist
    if pos == -1:
        if jPos != -1:
            return('Only J', jPos, curPh, foundC)
        else:
            return('No V/J found', -1, curPh, foundC)
    else:
        if jPos == -1:
            return('Only V',pos, curPh, foundC)
    return(aaSeq[pos:jPos], pos, curPh, foundC)

# Checks that the conserved amino acids remain
def isLegal(subAA):
    if (len(subAA)<4):
        return False
    if (subAA[0] == 'W'):
        if ((subAA[1] == 'G') & (subAA[3] == 'G')):
            return True
    if (subAA[0] == 'F'):
        if ((subAA[1] == 'G') | (subAA[3] == 'G')):
            return True
    if ((subAA[1] == 'G') & (subAA[3] == 'G')):
        return True
    if ((subAA[0] == 'G') & (subAA[2] == 'G')):
        return True
    if subAA.startswith('FS'):
        return True
    return False




def getNTseq(fullSeq):
    mod = len(fullSeq) % 3
    if mod != 0:
        fSeq = fullSeq[:-mod].translate()
    else:
        fSeq = fullSeq.translate()
    return fSeq

def findSeqAndLengthOfAA(aaSeq):
    fLen = 0
    fSeq = ''
    startArr = []
    stopArr = []
    startM = aaSeq.find('M')
    while startM != -1:
        startArr.append(startM)
        startM = aaSeq.find('M', startM + 1)
    stopPos = aaSeq.find('*')
    while stopPos != -1:
        stopArr.append(stopPos)
        stopPos = aaSeq.find('*', stopPos + 1)

    if ((len(startArr) == 0) | (len(stopArr) == 0)):
        return(fSeq,fLen)
    for stP in startArr:
        currStop = findStop(stP, stopArr)
        if currStop == -1:
            return (fSeq, fLen)
        else:
            currLen = currStop - stP
            if currLen >= fLen:
                fLen = currLen
                fSeq = aaSeq[stP:currStop]
    return fSeq


def findStop(stP, stopArr):
    for x in stopArr:
        if x > stP:
            return x
    return -1


def makeAADict(aaF):
    fDict = dict()
    f = open(aaF,'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        if lArr[0] in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Warning! %s appear twice in AA file\n' % lArr[0])
            sys.stderr.flush()
        fDict[lArr[0]] = lArr[1]
        l = f.readline()
    f.close()
    return fDict



def pickFinalIsoforms(fullTcrFileHeavy, fullTcrFileKappa,fullTcrFileLambda , output):
    pickFinalIsoformChain(fullTcrFileHeavy, output + '.heavy.full.BCRs.bestIso.fa', output + '.heavy.rsem.out.genes.results')
    pickFinalIsoformChain(fullTcrFileKappa, output + '.kappa.full.BCRs.bestIso.fa', output + '.kappa.rsem.out.genes.results')
    pickFinalIsoformChain(fullTcrFileLambda, output + '.lambda.full.BCRs.bestIso.fa', output + '.lambda.rsem.out.genes.results')



def pickFinalIsoformChain(fullTCRfa, newFasta, rsemF):
    if os.path.isfile(fullTCRfa):
        f = open(fullTCRfa, 'rU')
        outFa = open(newFasta, 'w')
        fastaDict = dict()
        byVJDict = dict()
        for record in SeqIO.parse(f,'fasta'):
            if record.id in fastaDict:
            #record.id = record.id + '_2'
                sys.stderr.write(str(datetime.datetime.now()) + 'error! same name for two fasta entries %s\n' % record.id)
                sys.stderr.flush()
            fastaDict[record.id] = record.seq
            onlyVJrec = str(record.id)
            idArr = onlyVJrec.strip('\n').split('.')
            vjStr = idArr[3] + '.' + idArr[4]
            if vjStr not in byVJDict:
                byVJDict[vjStr] = []
            byVJDict[vjStr].append(record.id)
        for vjStr in byVJDict:

            if len(byVJDict[vjStr]) == 1:
                cId = byVJDict[vjStr][0]
                cSeq = fastaDict[cId]
                newRec = SeqRecord(cSeq, id = cId, description = '')
                SeqIO.write(newRec,outFa,'fasta')
            else:
            #print vjStr
            #print byVJDict[vjStr]
                bestId = findBestC(byVJDict[vjStr], rsemF)
            #print "best: " + bestId
                bSeq = fastaDict[bestId]
                newRec = SeqRecord(bSeq, id = bestId, description = '')
                SeqIO.write(newRec,outFa,'fasta')
        outFa.close()
        f.close()

def findBestC(vjArr, rsemF):
    if (os.path.exists(rsemF)):
        f = open(rsemF, 'r')
        f.readline()
        l = f.readline()
        bestSeq = 'name'
        maxCount = 0.0
        while l != '':
            lArr = l.strip('\n').split('\t')
            if lArr[0] in vjArr:
                currCount = float(lArr[4])
                if currCount > maxCount:
                    bestSeq = lArr[0]
                    maxCount = currCount
            l = f.readline()
        f.close()
        if bestSeq == 'name':
            return vjArr[0]
        return bestSeq
    else:
        return vjArr[0]


def runRsem(outDir, rsem, bowtie2, fullTcrFileHeavy, fullTcrFileKappa,fullTcrFileLambda, output, samtools, secondRound):
    if samtools != '':
        if samtools[-1] != '/':
            samtools += '/'
    rsemIndDir = outDir + 'rsem_ind' + str(random.randint(1,100000000))
    if os.path.exists(rsemIndDir) == False:
        os.makedirs(rsemIndDir)
    if rsem != '':
        if rsem[-1] != '/':
            rsem += '/'
    if bowtie2 != '':
        if bowtie2[-1] != '/':
            bowtie2 += '/'
    runRSEMonOneFile(fullTcrFileHeavy, rsem, bowtie2, rsemIndDir, output, samtools, "heavy", secondRound)
    runRSEMonOneFile(fullTcrFileKappa, rsem, bowtie2, rsemIndDir, output, samtools, "kappa", secondRound)
    runRSEMonOneFile(fullTcrFileLambda, rsem, bowtie2, rsemIndDir, output, samtools, "lambda", secondRound)
    shutil.rmtree(rsemIndDir)



def runRSEMonOneFile(resF, rsem, bowtie2, rsemIndDir, output, samtools, chain, secondRound):
    if os.path.exists(resF):
        if bowtie2 != '':
            if not secondRound:
                with open(os.devnull, 'w') as devnull:
                    subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2 ,
                                     '-q', resF, rsemIndDir + '/VDJ.' + chain + '.seq'], stdout=devnull, stderr=devnull)
                    subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                                     '--bowtie2', '--bowtie2-path',bowtie2, '--bowtie2-mismatch-rate', '0.0' , '--paired-end',
                                     output + '.' + chain + '.R1.fa', output + '.' + chain + '.R2.fa',
                                     rsemIndDir + '/VDJ.' + chain + '.seq', output + '.' + chain + '.rsem.out'], stdout=devnull, stderr=devnull)
            else:
                with open(os.devnull, 'w') as devnull:
                    subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2', '--bowtie2-path', bowtie2 ,
                                     '-q', resF, rsemIndDir + '/VDJ.' + chain + '.seq'], stdout=devnull, stderr=devnull)
                    subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                                     '--bowtie2', '--bowtie2-path',bowtie2, '--bowtie2-mismatch-rate', '0.0' , '--paired-end',
                                     output + '.' + chain + '.R1.fa', output + '.' + chain + '.R2.fa',
                                     rsemIndDir + '/VDJ.' + chain + '.seq', output + '.afterCDRrec.' + chain + '.rsem.out'], stdout=devnull, stderr=devnull)

        else:
            if not secondRound:
                with open(os.devnull, 'w') as devnull:
                    subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                                     '-q', resF, rsemIndDir + '/VDJ.' + chain + '.seq'], stdout=devnull, stderr=devnull)
                    subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                                     '--bowtie2', '--bowtie2-mismatch-rate', '0.0', '--paired-end',
                                     output + '.' + chain + '.R1.fa', output + '.' + chain + '.R2.fa',
                                     rsemIndDir + '/VDJ.' + chain + '.seq', output + '.' + chain + '.rsem.out'], stdout=devnull, stderr=devnull)
            else:
                with open(os.devnull, 'w') as devnull:
                    subprocess.call([rsem + 'rsem-prepare-reference' , '--bowtie2',
                                     '-q', resF, rsemIndDir + '/VDJ.' + chain + '.seq'], stdout=devnull, stderr=devnull)
                    subprocess.call([rsem + 'rsem-calculate-expression', '--no-qualities', '-q',
                                     '--bowtie2', '--bowtie2-mismatch-rate', '0.0', '--paired-end',
                                     output + '.' + chain + '.R1.fa', output + '.' + chain + '.R2.fa',
                                     rsemIndDir + '/VDJ.' + chain + '.seq', output + '.afterCDRrec.' + chain + '.rsem.out'], stdout=devnull, stderr=devnull)
        if not secondRound:
            unsortedBam = output + '.' + chain + '.rsem.out.transcript.bam'
        else:
            unsortedBam = output + '.afterCDRrec.' + chain + '.rsem.out.transcript.bam'
        if not os.path.exists(unsortedBam):
            print("RSEM did not produce any transcript alignment files for " + chain + " chain, please check the -rsem parameter")
        else:
            if not secondRound:
                sortedBam = output + '.' + chain + '.rsem.out.transcript.sorted.bam'
            else:
                sortedBam = output + '.afterCDRrec.' + chain + '.rsem.out.transcript.sorted.bam'
            if not os.path.exists(sortedBam):
                print("sorting bam file")
                subprocess.call([samtools + 'samtools', 'sort','-o',sortedBam, unsortedBam])
                subprocess.call([samtools + 'samtools', 'index', sortedBam])
    else:
        sys.stdout.write(str(datetime.datetime.now()) + " Did not reconstruct any " + chain + " chains, not running RSEM on this chain\n")
        sys.stdout.flush()


def createTCRFullOutput(fastaDict, tcr, outName, bases, mapDict, cSeq, cName, cId, oneSide):
    tcrF = open(tcr, 'rU')
    found = False
    ffound = False
    recNameArr = []
    for tcrRecord in SeqIO.parse(tcrF, 'fasta'):
        addedC = False
        tcrSeq = str(tcrRecord.seq)
        if tcrSeq.find('NNNNN') == -1 :
            if ffound == False:
                ffound = True
                outF = open(outName, 'w')
            idArr = tcrRecord.id.split('.')
            vEns = idArr[0]
            jEns = idArr[1].split('(')[0]
            vSeq = fastaDict[vEns]
            jSeq = fastaDict[jEns]
            recNameArr = writeRecord(tcrRecord, tcrSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict,bases, cSeq, cId, cName, outF,fastaDict, recNameArr, False)
    tcrF.close()
    tcrF = open(tcr, 'rU')
    for tcrRecord in SeqIO.parse(tcrF, 'fasta'):
        addedC = False
        tcrSeq = str(tcrRecord.seq)
        if ((tcrSeq.find('NNNNN') != -1) & (oneSide)) :
            curSeq = tcrSeq.split('NNNNNN')[0]
            jSeg = findBestJforSeq(curSeq,fastaDict,mapDict)
            if jSeg != 'NA':
                if ffound == False:
                    ffound = True
                    outF = open(outName, 'w')
                idArr = tcrRecord.id.split('.')
                vEns = idArr[0]
                vSeq = fastaDict[vEns]
                for jEns in jSeg:
                    jSeq = fastaDict[jEns]
                    recNameArr = writeRecord(tcrRecord, curSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict,bases, cSeq, cId, cName, outF,fastaDict, recNameArr, True)
            curSeqArr = tcrSeq.split('NNNNNN')
            curSeq = curSeqArr[-1]
            vSeg = findBestVforSeq(curSeq,fastaDict,mapDict)
            #print "found V segment:\n"
            #for v in vSeg:
            #    print v + '\n'
            if vSeg != 'NA':
                if ffound == False:
                    ffound = True
                    outF = open(outName, 'w')
                idArr = tcrRecord.id.split('.')
                jEns = idArr[1].split('(')[0]
                jSeq = fastaDict[jEns]
                for vEns in vSeg:
                    vSeq = fastaDict[vEns]
                    recNameArr = writeRecord(tcrRecord, curSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict,bases, cSeq, cId, cName, outF,fastaDict, recNameArr, True)

    tcrF.close()
    if found == True:
        outF.close()

def writeRecord(tcrRecord, tcrSeq, addedC, vEns, jEns, vSeq, jSeq, mapDict, bases, cSeq, cId, cName, outF,fastaDict, recNameArr, oneSide):
    vSeqTrim = ''
    jSeqTrim = ''
    if bases == -10:
        bases = min(len(vSeq), len(jSeq))
    elif bases > len(jSeq):
        jSeq = jSeq + cSeq
        addedC = True
    found = False
    for i in reversed(range(20,bases)):
        juncStart = tcrSeq[:i]
        vInd = vSeq.find(juncStart)
        if (vInd != -1):
            found = True
            vSeqTrim = vSeq[:vInd]
            break
    if found == False:
        vSeqTrim = vSeq[:-bases]
    found = False
    for j in reversed(range(20,bases)):
        juncEnd = tcrSeq[-j:]
        jInd = jSeq.find(juncEnd)
        if (jInd != -1):
            found = True
            jSeqTrim = jSeq[jInd + j:]
            break
    if found == False:
        jSeqTrim = jSeq[bases:]
    # Add TRBC or TRAC
    cArr = []
    if (str(tcrRecord.id).find('IGH')!= -1):
        for ens in mapDict:
            if mapDict[ens] in ['IGHA','IGHE','IGHG2C','IGHG2B','IGHG2A','IGHG1','IGHG3','IGHD','IGHM']:
                cArr.append(ens)
    elif (str(tcrRecord.id).find('IGK')!= -1):
        for ens in mapDict:
            if mapDict[ens].find('IGKC') != -1:
                cArr.append(ens)
    else:
        for ens in mapDict:
            if mapDict[ens].find('IGLC') != -1:
                cArr.append(ens)
    if not addedC:
        for ens in cArr:
            if ens in fastaDict:
                cSeq = fastaDict[ens]
                newSeq = vSeqTrim + tcrSeq + jSeqTrim + cSeq
                newId = mapDict[vEns] + '.' + mapDict[jEns] + '.' + mapDict[ens] + '.' + vEns + '.' + jEns + '.' + ens
                # Make sure we're not writing to BCRs with the same V/J/C combination, one from a full reconstruction and
                # one from a oneSided reconstruction
                if not ((newId in recNameArr) & (oneSide) ):
                    while newId in recNameArr:
                        newId += '_2'
                    recNameArr.append(newId)
                    record = SeqRecord(Seq(newSeq,IUPAC.ambiguous_dna), id = newId, description = '')
                    SeqIO.write(record,outF,'fasta')
    else:
        newSeq = vSeqTrim + tcrSeq + jSeqTrim
        newId = mapDict[vEns] + '.' + mapDict[jEns] + '.' + cName + '.' + vEns + '.' + jEns + '.' + cId
        if not ((newId in recNameArr) & (oneSide) ):
            while newId in recNameArr:
                newId += '_2'
            recNameArr.append(newId)
            record = SeqRecord(Seq(newSeq,IUPAC.ambiguous_dna), id = newId, description = '')
            SeqIO.write(record,outF,'fasta')
    return recNameArr



def findBestJforSeq(curSeq,fastaDict,idNameDict):
    jArrOld = findJsPerLen(curSeq, fastaDict, idNameDict,20)
    if len(jArrOld) == 0:
        return 'NA'
    for x in range(21,len(curSeq)):
        newArr = findJsPerLen(curSeq, fastaDict, idNameDict,x)
        if len(newArr) == 0:
            return jArrOld
        else:
            jArrOld = newArr
    print('Found a full J segment as the V/J junction, ignoring this reconstruction')
    return 'NA'

def findBestVforSeq(curSeq,fastaDict,idNameDict):
    vArrOld = findVsPerLen(curSeq, fastaDict, idNameDict,18)
    if len(vArrOld) == 0:
        return 'NA'
    for x in range(19,len(curSeq)):
        newArr = findVsPerLen(curSeq, fastaDict, idNameDict,x)
        if len(newArr) == 0:
            return vArrOld
        else:
            vArrOld = newArr
    print('Found a full V segment as the V/J junction, ignoring this reconstruction')
    return 'NA'

def findVsPerLen(curSeq, fastaDict, idNameDict,trim):
    fArr = []
    for seq in fastaDict:
        if idNameDict[seq].find('V') != -1:
            vSeq = fastaDict[seq]
            lenV = len(vSeq)
            for i in reversed(range(0,lenV+1)):
                if ((i - trim) >= 0):
                    #change this:
                    trimV = vSeq[i-trim:i]
                    if curSeq.find(trimV) != -1:
                        if seq not in fArr:
                            fArr.append(seq)
                            break
    return fArr





def findJsPerLen(curSeq, fastaDict, idNameDict,trim):
    fArr = []
    for seq in fastaDict:
        if idNameDict[seq].find('J') != -1:
            jSeq = fastaDict[seq]
            lenJ = len(jSeq)
            for i in range(0,lenJ):
                if ((i + trim) <= lenJ):
                    trimJ = jSeq[i:i+trim]
                    if curSeq.find(trimJ) != -1:
                        if seq not in fArr:
                            fArr.append(seq)
                            break
    return fArr




def analyzeChain(fastaDict, vdjDict, output, bam, unmapped, idNameDict, bases, chain, strand, NolowQ, top, byExp, readOverlap, downsample, organism, Hminus, Kminus, Lminus):
    junctionSegs = makeJunctionFile(bam, chain, output, bases, vdjDict, fastaDict, idNameDict, top, byExp, readOverlap, organism)
    unDict = writeReadsFile(bam, unmapped, junctionSegs, output, vdjDict, chain, strand, NolowQ, downsample, organism, Hminus, Kminus, Lminus)
    return unDict

def getCInfo(bedEntry, idNameDict, fastaDict):
    bedArr = bedEntry.strip('\n').split('\t')
    cId = bedArr[3]
    cName = idNameDict[cId]
    cSeq = fastaDict[cId]
    return(cSeq,cName,cId)


def makeJunctionFile(bam, chain, output, bases, vdjDict, fastaDict, idNameDict, top, byExp, readOverlap, organism):
    mappedFile = pysam.AlignmentFile(bam,"rb")
    if chain == 'H':
        vdjChainDict = vdjDict['heavy']
        outName = output + '.heavy.junctions.txt'
    elif chain == 'K':
        vdjChainDict = vdjDict['kappa']
        outName = output + '.kappa.junctions.txt'
    else:
        vdjChainDict = vdjDict['lambda']
        outName = output + '.lambda.junctions.txt'
    jSegs = vdjChainDict['J']
    vSegs = vdjChainDict['V']
    (cSeq, cId, cName) = getCInfo(vdjChainDict['C'][0], idNameDict, fastaDict)
    vjSegs = []
    for x in jSegs:
        vjSegs.append(x)
    for y in vSegs:
        vjSegs.append(y)
    vjReads = dict()
    (vjReads, vjCounts) = loadReadsToDict(vjSegs, mappedFile, vjReads, readOverlap)
    junctionSegs = writeJunctions(vjReads,outName, bases, fastaDict, idNameDict, cSeq, top, vjCounts, byExp)
    if len(junctionSegs) == 0:
        sys.stdout.write(str(datetime.datetime.now()) + ' Did not find any V-J reads, searching for V-C and J-C reads:\n')
        sys.stdout.flush()
        cReads = dict()
        (cReads, cCountsDict) = loadReadsToDict(vdjChainDict['C'], mappedFile, cReads, readOverlap)
        junctionSegs = writeJunctionsWithC(vjReads,outName, bases, fastaDict, idNameDict, cReads, cSeq)
    mappedFile.close()
    return junctionSegs

# Similar to "writeJunctionsWithC, only that instead of looking for V-J paired-reads, it looks for
# V-C and J-C paired-reads
# INPUT:
#       vjReads - reads dict of the V and J segments created by loadReadsToDict
#       outName - output name for junction file
#       bases - number of bases to take from V and J for the junction
#       fastaDict
#       idNameDict
#       cReads - reads dict of the C segments created by loadReadsToDict
# OUTPUT:
#        fArr - the V and J segments for which we found a junction for
def writeJunctionsWithC(vjReads,outName, bases, fastaDict, idNameDict, cReads, cSeqFa):
    out = open(outName,'w')
    fArr = []
    vArr = []
    jArr = []
    for seg in vjReads:
        if len(vjReads[seg]) > 0 :
            for cSeg in cReads:
                if (len(cReads[cSeg]) > 0) :
                    if (len([val for val in vjReads[seg]['first'] if val in cReads[cSeg]['second']]) > 0) |\
                                (len([val for val in vjReads[seg]['second'] if val in cReads[cSeg]['first']]) > 0) :
                        if idNameDict[seg].find('J') != -1 :
                            if seg not in jArr:
                                jArr.append(seg)
                        elif idNameDict[seg].find('V') != -1 :
                            if seg not in vArr:
                                vArr.append(seg)
                        fArr.append(seg)
    for vSeg in vArr:
        for jSeg in jArr:
            vSeqFa = fastaDict[vSeg]
            jSeqFa = fastaDict[jSeg]
            if len(jSeqFa) < bases:
                jSeqFa = jSeqFa + cSeqFa
            lenSeg = min(len(vSeqFa),len(jSeqFa))
            if bases != -10:
                if lenSeg < bases:
                    if len(vSeqFa) < len(jSeqFa):
                        lenSeg = len(vSeqFa)
                    else:
                        lenSeg = len(jSeqFa)
                    sys.stdout.write(str(datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V or J segment, taking the length' \
                                        'of the V/J segment instead, which is: ' + str(lenSeg) + '\n')
                    sys.stdout.flush()
                else:
                    lenSeg = bases
            jTrim = jSeqFa[:lenSeg]
            vTrim = vSeqFa[-1*lenSeg:]
            junc = vTrim + jTrim
            recordName = vSeg + '.' + jSeg + '(' + idNameDict[vSeg] + '-' + idNameDict[jSeg] + ')'
            record = SeqRecord(Seq(junc,IUPAC.ambiguous_dna), id = recordName, description = '')
            SeqIO.write(record,out,'fasta')
    out.close()
    return fArr


# Load all the reads from a list of segments into a dictionary.
# INPUT: segsDict: A dict, where the key is a segment, and the value is an array of bed entries of this segment
#        mappedFile: Bam file of the mapped reaeds
#       readDict: A dictionary. The keys are the segment name, the values is a dictionary 'first':[] and 'second':[]
#                 where 'first' array holds the query name of R1's that overlap this segment, and 'second' holds the
#                 query name of R2's that overlap the segment.
def loadReadsToDict(segsDict, mappedFile, readDict, readOverlap):
    countDict = dict()
    for seg in segsDict:
        lArr = seg.strip('\n').split('\t')
        segName = lArr[3]
        readDict[segName] = {'first':[],'second':[]}
        lArr = seg.strip('\n').split('\t')
        chr = lArr[0]
        start = int(lArr[1])
        end = int(lArr[2])
        readsIter = mappedFile.fetch(chr, start-1, end+1)
        readCounter = 0
        for read in readsIter:
            #make sure the read overlap at least "readOverlap" bases of the segment.
            overlap = read.get_overlap(start-1,end+1)
            if (end-start) < readOverlap:
                readOverlap = end-start-15
            if overlap >= readOverlap:
                currName = read.query_name
                if read.is_read1:
                    if currName not in readDict[segName]['first']:
                        readDict[segName]['first'].append(currName)
                        readCounter += 1
                elif read.is_read2:
                    if currName not in readDict[segName]['second']:
                        readDict[segName]['second'].append(currName)
                        readCounter += 1
        countDict[segName] = readCounter
    return (readDict, countDict)


def writeReadsFile(bam, unmapped, junctionSegs, output, vdjDict, chain, strand, NolowQ, downsample, organism, Hminus, Kminus, Lminus):
    if chain == 'H':
        vdjChainDict = vdjDict['heavy']
        outReads = output + '.heavy.mapped.and.unmapped.fa'
        pairedReads1 = output + '.heavy.R1.fa'
        pairedReads2 = output + '.heavy.R2.fa'
    elif chain == 'K':
        vdjChainDict = vdjDict['kappa']
        outReads = output + '.kappa.mapped.and.unmapped.fa'
        pairedReads1 = output + '.kappa.R1.fa'
        pairedReads2 = output + '.kappa.R2.fa'
    else:
        vdjChainDict = vdjDict['lambda']
        outReads = output + '.lambda.mapped.and.unmapped.fa'
        pairedReads1 = output + '.lambda.R1.fa'
        pairedReads2 = output + '.lambda.R2.fa'
    out = open(outReads, 'w')
    constDict = vdjChainDict['C']
    # This dict for unmapped reads has reads that should be rev.comp in the revcomp arr, otherwise in id.
    # For every read, the value is a tuple - the first value is first/second, to make sure there are no errors.
    # The second value is id/revcomp, to see what sequence should be written.
    # Note: This classification is about what should happen to the unmapped reads, not how their paired mapped
    # reads were read.
    unmappedDict = dict()
    seqDict = dict()
    oriDict = dict()
    alignedDict = dict()
    lowQDict = dict()
    primDict = dict()
    for seg in constDict:
        (unmappedDict, alignedDict, seqDict, lowQDict, oriDict, primDict) = addReadsToDict(unmappedDict, seg, bam,  False, alignedDict, seqDict, strand, 'C', lowQDict, chain, oriDict, primDict, organism)
    vSegs = vdjChainDict['V']
    for vSeg in vSegs:
        vSegName = vSeg.strip('\n').split('\t')[3]
        if vSegName in junctionSegs:
            (unmappedDict, alignedDict, seqDict, lowQDict, oriDict, primDict) = addReadsToDict(unmappedDict, vSeg, bam,  True, alignedDict, seqDict, strand, 'V', lowQDict, chain, oriDict, primDict, organism)
    jSegs = vdjChainDict['J']
    for jSeg in jSegs:
        jSegName = jSeg.strip('\n').split('\t')[3]
        if jSegName in junctionSegs:
            (unmappedDict, alignedDict, seqDict, lowQDict, oriDict, primDict) = addReadsToDict(unmappedDict, jSeg, bam,  True, alignedDict, seqDict, strand, 'J', lowQDict, chain, oriDict, primDict, organism)
    unDict = dict()
    fullAlignDict = alignedDict
    if downsample:
        if len(alignedDict) > 5000:
            curLen = len(alignedDict)
            sys.stdout.write(str(datetime.datetime.now()) + ' Found %s reads mapped to V/J segments, downsampling to 5000 for reconstruction\n' % str(curLen))
            sys.stdout.flush()
            alignedDict = {k: v for k, v in random.sample(alignedDict.items(),5000)}
    for read in alignedDict:
        record = SeqRecord(Seq(alignedDict[read], IUPAC.ambiguous_dna), id = read, description = '')
        SeqIO.write(record,out,'fasta')
    (seqDict,unDict) = writeUnmappedReads(unmappedDict, out, unmapped, seqDict, unDict, fullAlignDict, lowQDict, NolowQ, chain, downsample, organism, Hminus, Kminus, Lminus)
    seqDict = addMappedPairsToSeqDict(seqDict, bam, out, NolowQ, fullAlignDict, oriDict)
    writeSeqDict(seqDict, pairedReads1, pairedReads2)
    out.close()
    return unDict

def addMappedPairsToSeqDict(seqDict, bam, out, NolowQ, alignedDict, oriDict):
    firstDict = dict()
    secondDict = dict()
    for name in seqDict:
        # the value of firstDict[name] or secondDict[name] is the orientation of the mapped pair. When writing the new
        # read we should write the opposite of the value stated in the dict
        if ((seqDict[name][1] == '1') & (seqDict[name][0] == '0')):
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! empty record insdie seqDict\n')
            sys.stderr.flush()
        if seqDict[name][0] == '0':
            firstDict[name] = oriDict[name][1]
            if oriDict[name][1] == 'NA':
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! oriDict is NA when seqDict is not\n')
                sys.stderr.flush()
            if oriDict[name][0] != 'NA':
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconrd have no seq in seqDict but have orientation in oriDict\n')
                sys.stderr.flush()
        if seqDict[name][1] == '1':
            secondDict[name] = oriDict[name][0]
            if oriDict[name][0] == 'NA':
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! oriDict is NA when seqDict is not\n')
                sys.stderr.flush()
            if oriDict[name][1] != 'NA':
                sys.stderr.write(str(datetime.datetime.now()) + ' Error! reconrd have no seq in seqDict but have orientation in oriDict\n')
                sys.stderr.flush()
    f = pysam.AlignmentFile(bam,"rb")
    readsIter = f.fetch()
    for read in readsIter:
        name = read.query_name
        pos = -1
        mateOri = 'NA'
        if ((name in firstDict) & (read.is_read1)):
            mateOri = firstDict[name]
            pos = 0
            check = '0'
        elif ((name in secondDict) & (read.is_read2)):
            mateOri = secondDict[name]
            pos = 1
            check = '1'
        if pos != -1:
            qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
            if read.is_read1:
                currName = name + '\\1'
            else:
                currName = name + '\\2'
            if not NolowQ:
                if mateOri == 'id':
                    if read.is_reverse:
                        rSeq = qSeq
                    else:
                        rSeq = qSeq.reverse_complement()
                elif mateOri == 'rev':
                    if read.is_reverse:
                        rSeq = qSeq.reverse_complement()
                    else:
                        rSeq = qSeq
                else:
                    sys.stderr.write(str(datetime.datetime.now()) + ' Error! mate orientation is not rev or id\n')
                    sys.stderr.flush()
                if currName not in alignedDict:
                    alignedDict[currName] = str(rSeq)
                    record = SeqRecord(rSeq, id = currName, description = '')
                    SeqIO.write(record,out,'fasta')
                else:
                    if str(alignedDict[currName]) != str(rSeq):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Error! read %s has two different sequences in alignedDict\n' % currName)
                        sys.stderr.write(str(datetime.datetime.now()) + ' First: %s\n' % alignedDict[currName])
                        sys.stderr.write(str(datetime.datetime.now()) + ' Second: %s\n' % str(rSeq))
                        sys.stderr.flush()
            if read.is_reverse:
                qSeq = qSeq.reverse_complement()
            if seqDict[name][pos] != check:
                if str(qSeq) != str(seqDict[name][pos]):
                    sys.stderr.write(str(datetime.datetime.now()) + ' Error! read %s has two different mapped sequences not in the V/J region\n' % name)
                    sys.stderr.flush()
            seqDict[name][pos] = qSeq
    f.close()
    return seqDict


def writeSeqDict(seqDict, r1, r2):
    r1f = open(r1,'w')
    r2f = open(r2,'w')
    for seq in seqDict:
        if ((seqDict[seq][0] != '0') & (seqDict[seq][1]!= '1')):
            seq1 = seq
            seq2 = seq
            rec1 = SeqRecord(seqDict[seq][0], id = seq1, description = '')
            rec2 = SeqRecord(seqDict[seq][1], id = seq2, description = '')
            SeqIO.write(rec1,r1f,'fasta')
            SeqIO.write(rec2,r2f,'fasta')
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' The read %s has only one mate found, ignoring it\n' % seq)
            sys.stderr.flush()
    r1f.close()
    r2f.close()

def writeUnmappedReads(unmappedDict, out, unmapped, seqDict, unDict, alignedDict, lowQDict, NolowQ, chain, downsample, organism, Hminus, Kminus, Lminus):
    count = 0
    if downsample:
        if len(unmappedDict) > 5000:
            curLen = len(unmappedDict)
            sys.stdout.write(str(datetime.datetime.now()) + ' Found %s unmapped reads, downsampling to 5000 for reconstruction\n' % str(curLen))
            sys.stdout.flush()
            unmappedDictToWrite = {k: v for k, v in random.sample(unmappedDict.items(),5000)}
        else:
            unmappedDictToWrite = unmappedDict
    else:
        unmappedDictToWrite = unmappedDict
    f = pysam.AlignmentFile(unmapped,"rb")
    readsIter = f.fetch(until_eof = True)
    for read in readsIter:
        name = read.query_name
        if name in unmappedDict:
            cName = name
            unDictName = name
            (strand , ori) = unmappedDict[name]
            if (((strand == 'first') & (read.is_read2) ) | ((strand == 'second') & (read.is_read1))):
                sys.stderr.write(str(datetime.datetime.now()) + ' Warning! unmapped read %s is inconsistent regarding first/second read\n' % cName)
                sys.stderr.flush()
            else:
                if strand == 'first' :
                    name += '\\1'
                    unDictName += '_1'
                else:
                    name += '\\2'
                    unDictName += '_2'
                if unDictName in unDict:
                    sys.stderr.write(str(datetime.datetime.now()) + ' Warning! unmapped read %s appear twice in unmapped bam file\n' % cName)
                    sys.stderr.flush()
                unDict[unDictName] = '1'
                qSeq = Seq(read.query_sequence, IUPAC.ambiguous_dna)
                if ori == 'rev':
                    qSeq = qSeq.reverse_complement()
                if organism in ['mm10','mm10_ncbi','hg38']:
                    if ((chain == 'H') | ((chain == 'L') & (organism in ['mm10','mm10_ncbi']))):
                        qSeq = qSeq.reverse_complement()
                else:
                    if (((chain == 'H') & (Hminus)) | ((chain == 'L') & (Lminus)) | ((chain == 'K') & (Kminus))):
                        qSeq = qSeq.reverse_complement()
                if name in alignedDict:
                    if alignedDict[name] != str(qSeq):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Warning! unmapped read %s appear twice in alignedDict with differnet seqs\n' % name)
                        sys.stderr.flush()
                else:
                    if ((not NolowQ) | ((name not in lowQDict) & (NolowQ))):
                        alignedDict[name] = str(qSeq)
                        if cName in unmappedDictToWrite:
                            count += 1
                            record = SeqRecord(qSeq, id = name, description = '')
                            SeqIO.write(record,out,'fasta')
                if ori == 'rev':
                    qSeq = qSeq.reverse_complement()
                if cName not in seqDict:
                    sys.stderr.write(str(datetime.datetime.now()) + ' Warning! unmapped read is in unmappedDict but not in seqDict %s\n' % read.query_name)
                    sys.stderr.flush()
                else:
                    if strand == 'first':
                        seqDict[cName][0] = qSeq
                    else:
                        seqDict[cName][1] = qSeq
    f.close()
    print("Added " + str(count) + " unmapped reads to the reads file")
    return (seqDict, unDict)

# Aligned dict - all the reads (with _1/_2) that were already written to the mapped.unmapped.fa file
def addReadsToDict(unmappedDict, segBed, bam, mappedRead, alignedDict, seqDict, strand, segType, lowQDict, chain, oriDict, primDict, organism):
    bedArr = segBed.strip('\n').split('\t')
    chr = bedArr[0]
    start = int(bedArr[1])
    end = int(bedArr[2])
    mappedFile = pysam.AlignmentFile(bam,"rb")
    readsIter = mappedFile.fetch(chr, start-1, end+1)
    for read in readsIter:
        currName = read.query_name
        if currName not in seqDict:
            seqDict[currName] = ['0','1']
            oriDict[currName] = ['NA','NA']
            primDict[currName] = ['NA','NA']
        if read.is_read1:
            pairRead = 'second'
            readName = currName + '\\1'
            pairName = currName + '\\2'
            seqPos = 0
        elif read.is_read2:
            pairRead = 'first'
            readName = currName + '\\2'
            pairName = currName + '\\1'
            seqPos = 1
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read is not read1 and not read2\n')
            sys.stderr.flush()
        currSeq = Seq(read.query_sequence,IUPAC.ambiguous_dna)
        if read.is_reverse:
            pairOr = 'id'
            oriDict[currName][seqPos] = 'rev'
            readStrand = 'minus'
            currSeq = currSeq.reverse_complement()
        else:
            readStrand = 'plus'
            pairOr = 'rev'
            oriDict[currName][seqPos] = 'id'
        seqDict[currName][seqPos] = currSeq
        if read.mate_is_unmapped:
            takePair = toTakePair(segType, strand, readStrand, chain, organism)
            if takePair == False:
                lowQDict[pairName] = '1'
            #takePair = True
            #if takePair:
            if currName in unmappedDict:
                if unmappedDict[currName] != (pairRead, pairOr):
                    if mappedRead:
                        if ((primDict[currName][seqPos] == 'primary') & (not read.is_secondary)):
                            #sys.stderr.write(str(datetime.datetime.now()) + ' PairRead: ' + pairRead + ' pair or: ' + pairOr + '\n' )
                            #sys.stderr.write(str(datetime.datetime.now()) + ' unmappedDict[currName]: ' + unmappedDict[currName][0] + ' ' + unmappedDict[currName][1] + '\n' )
                            sys.stderr.write(str(datetime.datetime.now()) + ' Warning! Read %s has more than one unmppaed mate with differnet strand/mate, both are primary alignments\n' % currName)
                            sys.stderr.flush()
                        elif ((primDict[currName][seqPos] == 'NA') & (not read.is_secondary)):
                            unmappedDict[currName] = (pairRead, pairOr)
                    else:
                        sys.stderr.write(str(datetime.datetime.now()) + ' Warning! Read %s has more than one unmppaed mate with differnet strand/mate, and mapping to C segment\n' % currName)
                        sys.stderr.flush()
            else:
                unmappedDict[currName] = (pairRead, pairOr)

        if mappedRead == True :
            if readName in alignedDict:
                if (alignedDict[readName] != read.query_sequence):
                    if ((primDict[currName][seqPos] == 'primary') & (not read.is_secondary)):
                        sys.stderr.write(str(datetime.datetime.now()) + ' Warning! Read %s has two instances of primary alignments but different seuqences\n' % read.query_name)
                        sys.stderr.flush()
                    else:
                        if ((primDict[currName][seqPos] == 'NA') & (not read.is_secondary)):
                            alignedDict[readName] = read.query_sequence
                            primDict[currName][seqPos] = 'primary'
            else:
                alignedDict[readName] = read.query_sequence
                if not read.is_secondary:
                    primDict[currName][seqPos] = 'primary'
    mappedFile.close()
    return(unmappedDict, alignedDict, seqDict, lowQDict, oriDict, primDict)


### For minus, add an unmapped pair if the current mate is: 1. V and Plus 2. J/C and minus
### For plus, exactly the opposite
### For the Heavy chain, do exactly the oppsite because its not the other strand

def toTakePair(segType, strand, readStrand, chain, organism):
    #TODO: testing this!!
    #chain = 'H'
    if strand == 'none':
        return True
    if ((readStrand != 'minus') & (readStrand != 'plus')):
        sys.stderr.write(str(datetime.datetime.now()) + ' Error! Read strand should be plus or minus only\n')
        sys.stderr.flush()
        return True
    if ((chain != 'H') & (not ((chain == 'L') & (organism in ['mm10','mm10_ncbi'])))):
        cond = True
    else:
        cond = False
    if ((segType == 'C') | (segType == 'J')):
        if strand == 'minus':
            if readStrand == 'minus':
                if cond:
                    return True
                else:
                    return False
            else:
                if cond:
                    return False
                else:
                    return True
        else:
            if readStrand == 'minus':
                if cond:
                    return False
                else:
                    return True
            else:
                if cond:
                    return True
                else:
                    return False
    else:
        if strand == 'minus':
            if readStrand == 'minus':
                if cond:
                    return False
                else:
                    return True
            else:
                if cond:
                    return True
                else:
                    return False
        else:
            if readStrand == 'minus':
                if cond:
                    return True
                else:
                    return False
            else:
                if cond:
                    return False
                else:
                    return True
    return True


def writeJunctions(vjReads,outName, bases, fastaDict, idNameDict, cSeq, top, vjCountsDict, byExp):
    if bases > 50:
        sys.stdout.write(str(datetime.datetime.now()) + ' Warning! Bases parameter might be larger than the length of the J segments,'
                                                        'appending the C segment to the J segment if needed\n')
        sys.stdout.flush()
    out = open(outName,'w')
    fArr = []
    pairCountDict = dict()
    seqsDict = dict()
    for seg in vjReads:
        if idNameDict[seg].find('J') != -1 :
            if len(vjReads[seg]) > 0 :
                for sSeg in vjReads:
                    if ((idNameDict[sSeg].find('V') != -1) & (len(vjReads[sSeg]) > 0)) :
                        if (len([val for val in vjReads[seg]['first'] if val in vjReads[sSeg]['second']]) > 0) |\
                                (len([val for val in vjReads[seg]['second'] if val in vjReads[sSeg]['first']]) > 0) :
                            if seg not in fArr:
                                fArr.append(seg)
                            if sSeg not in fArr:
                                fArr.append(sSeg)
                            vSeq = fastaDict[sSeg]
                            jSeq = fastaDict[seg]
                            lenSeg = min(len(vSeq),len(jSeq))
                            if bases != -10:
                                if lenSeg < bases:
                                    if bases > len(vSeq):
                                        sys.stdout.write(str(datetime.datetime.now()) + ' Bases parameter is bigger than the length of the V segment, taking the length' \
                                              'of the V segment instead, which is: ' + str(len(vSeq)) + '\n')
                                        sys.stdout.flush()
                                        lenSeg = len(vSeq)
                                        if lenSeg > len(jSeq):
                                            jSeq = jSeq + cSeq
                                    else:
                                        jSeq = jSeq + cSeq
                                        lenSeg = bases
                                else:
                                    lenSeg = bases
                            jTrim = jSeq[:lenSeg]
                            vTrim = vSeq[-1*lenSeg:]
                            junc = vTrim + jTrim
                            recordName = sSeg + '.' + seg + '(' + idNameDict[sSeg] + '-' + idNameDict[seg] + ')'
                            record = SeqRecord(Seq(junc,IUPAC.ambiguous_dna), id = recordName, description = '')
                            curCont = vjCountsDict[seg] + vjCountsDict[sSeg]
                            pairCountDict[str(record.seq)] = curCont
                            seqsDict[str(record.seq)] = record
    sorted_pairs = sorted(pairCountDict.items(), key=operator.itemgetter(1), reverse=True)
    if ((top == -1) | (top > len(sorted_pairs))):
        for rec in seqsDict:
            SeqIO.write(seqsDict[rec] ,out,'fasta')
    else:
        if not byExp:
            for i in range(0,top):
                SeqIO.write(seqsDict[sorted_pairs[i][0]],out,'fasta')
        else:
            wrote = 1
            SeqIO.write(seqsDict[sorted_pairs[0][0]],out,'fasta')
            curCount = sorted_pairs[0][1]
            wroteSecond = False
            for i in range(1,len(sorted_pairs)):
                if sorted_pairs[i][1] == curCount:
                    if not wroteSecond:
                        wroteSecond = True
                        SeqIO.write(seqsDict[sorted_pairs[i][0]],out,'fasta')
                        wrote += 1
                else:
                    curCount = sorted_pairs[i][1]
                    wroteSecond = False
                    SeqIO.write(seqsDict[sorted_pairs[i][0]],out,'fasta')
                    wrote += 1
                if wrote == top:
                    break

    out.close()
    return fArr





# Create a dict {'Alpha':{'C':[bed],'V':[bed],'J':[bed]}, 'Beta':{'C':[],'V':[],'J':[]}}
def makeVDJBedDict(bed,idNameDict):
    fDict = {'heavy':{'C':[],'V':[],'J':[]}, 'kappa':{'C':[],'V':[],'J':[]}, 'lambda':{'C':[],'V':[],'J':[]}}
    f = open(bed, 'r')
    l = f.readline()
    while l != '':
        lArr = l.strip('\n').split('\t')
        gID = lArr[3]
        gName = idNameDict[gID]
        chain = ''
        if (gName.startswith('IGH')):
            chain = 'heavy'
        elif (gName.startswith('IGK')):
            chain = 'kappa'
        elif (gName.startswith('IGL')):
            chain = 'lambda'
        else:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! found gene id that is not lambda/kappa/heavy: ' + gName)
            sys.stderr.flush()
        if gName.find('V') != -1:
            fDict[chain]['V'].append(l)
        elif gName.find('J') != -1:
            fDict[chain]['J'].append(l)
        else:
            if gName.find('IGLCOR') == -1:
                fDict[chain]['C'].append(l)
        l = f.readline()
    f.close()
    return fDict





# Creates a dictionary of ENSEMBL ID -> fasta sequence
def makeFastaDict(fasta):
    inF = open(fasta,'rU')
    fastaDict = dict()
    for record in SeqIO.parse(inF, 'fasta'):
        fastaDict[record.id] = str(record.seq)
    inF.close()
    return fastaDict


# Creates a dictionary of ENSEMBL ID -> Gene name
def makeIdNameDict(mapping):
    f = open(mapping, 'r')
    fDict = dict()
    linesArr = f.read().split('\n')
    linesArr = filter(None, linesArr)
    f.close()
    for line in linesArr:
        lineArr = line.split('\t')
        id = lineArr[0]
        name = lineArr[1]
        ind = name.find('Gene:')
        if ind != -1:
            name = name[ind+len('Gene:'):]
        if id in fDict:
            sys.stderr.write(str(datetime.datetime.now()) + ' Error! %s appear twice in mapping file\n' % id)
            sys.stderr.flush()
        fDict[id] = name
    return fDict


def checkParameters(strand, path, sumF, Hr, Fr):
    if strand.lower() not in ['none','minus','plus']:
        sys.exit("-strand should be one of: none, minus, plus")
    if ((Hr < 0) | (Hr > 1)):
        sys.exit("Mutation rate for complementary determining regions should be between 0 and 1")
    if ((Fr < 0) | (Fr > 1)):
        sys.exit("Mutation rate for framework regions should be between 0 and 1")
    if path == '':
        sys.exit("Must include the -path parameter")
    if sumF == '':
        sys.exit("Must include the -sumF parameter")
    if not os.path.isdir(path):
        sys.exit("%s path does not exists. Please check your -path parameter and run again" % path)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-genome','-g','-G', help='Alignment genome. Currently supported: mm10, mm10_ncbi and hg38. Other genomes require the user to prepare their own annotation files', required=True)
    parser.add_argument('-NoLowQ', help='add if you want to remove \"low quality\" reads as input to the reconstruction '
                                                        'algorithm', action='store_true')
    parser.add_argument('-oneSide', help='add if you want to observe reconstrctuion only from the V side', action='store_true')
    parser.add_argument('-path','-p','-P', help='The path for the data directory. Assumes that every subdirectory'
                                                    'is a single cell', default='')
    parser.add_argument('-sumF', help='prefix for summary outputs', default='')
    parser.add_argument('-bowtie2','-bw','-BW', help='Path to bowtie2. If not used assumes that bowtie2 is in the'
                                                 'default path', default = '')
    parser.add_argument('-rsem','-RSEM', help='Path to rsem. If not used assumes that rsem is in the'
                                                'default path', default = '')
    parser.add_argument('-strand', help='Strand of the right most read in genomic coordinates. Options are: [minus, plus, '
                                        'none]. Defualt is minus', default = 'minus')
    parser.add_argument('-output','-out','-o','-O', help='output prefix, relative to /path/singleCellFolder', required=True)
    parser.add_argument('-bam', help='Input bam alignment file, relative to /path/singleCellFolder/ if working on multiple files', default = './tophat_output/picard_output/sorted.bam')
    parser.add_argument('-unmapped','-u','-U', help='bam file of the unmapped reads, relative to /path/singleCellFolder/', default = './tophat_output/unmapped.bam')
    parser.add_argument('-bases','-b','-B', help='Number of bases to take from each V and J segments, default is min(len(V), len(J) ', type=int, default=-10)
    parser.add_argument('-iterations','-iter','-i','-I', help='Number of iterations for the reconstruction'
                                                              'algorithm, default is 20', type=int, default=20)
    parser.add_argument('-samtools', help='Path to samtools. If not used assumes that samtools is in the default path', default = '')
    parser.add_argument('-score','-sc','-SC', help='Alignment score threshold. Default is 15', type=int, default=15)
    parser.add_argument('-top','-t','-T', help='Take only the top x combination of V and J, based on the sum '
                                               'number of reads that map to both. Default is to take all', type=int, default=-1)
    parser.add_argument('-readOverlap','-ro','-readoverlap', help='Add a read to list of mapped reads only if it maps at least X bases'
                                               'to the V/J/C segment. Default is 1', type=int, default=1)
    parser.add_argument('-byExp', help='if using the Top option, add this tag if you want to take only two chains from each'\
                                                        'read count, until top is reached', action='store_true')
    parser.add_argument('-downsample', help='Add this flag in case of a very large files (>100,000 mapped V/J/C reads to subsample the data'\
                                                        'with this option, the total number of mapped reads output by BRAPeS is lower than the true measure.', action='store_true')
    parser.add_argument('-skipHVR', help='Add this flag if you want to skip reconstruction of CDR1, CDR2 and the framework regions', action='store_true')
    parser.add_argument('-Hr','-max_h_mutation_rate', help='Maximum allowed mutation rate for complementary determining regions (CDR1 and CDR2). '
                                               'Default is 0.35.', type=float, default=0.35)
    parser.add_argument('-Fr','-max_f_mutation_rate', help='Maximum allowed mutation rate for framework regions (FR1, FR2, FR3 and FR4). '
                                               'Default is 0.2.', type=float, default=0.2)
    parser.add_argument('-overlap','-ol','-OL', help='Number of minimum bases that overlaps V and J ends,'
                                                              'default is 10', type=int, default=10)
    parser.add_argument('-Hminus', help='Only when running BRAPeS on user defined genomes. Add this flag if the annotation of the V/J segmenets'\
                                                        'of the heavy chain are on the negative strand', action='store_true')
    parser.add_argument('-Kminus', help='Only when running BRAPeS on user defined genomes. Add this flag if the annotation of the V/J segmenets'\
                                                        'of the kappa chain are on the negative strand', action='store_true')
    parser.add_argument('-Lminus', help='Only when running BRAPeS on user defined genomes. Add this flag if the annotation of the V/J segmenets'\
                                                        'of the lambda chain are on the negative strand', action='store_true')

    args = parser.parse_args()
    runBCRpipe(args.genome, args.output, args.bam, args.unmapped, args.bases, args.strand,
                args.iterations,args.score, args.overlap, args.rsem, args.bowtie2,
                  args.path, args.sumF, args.NoLowQ, args.samtools, args.top, args.byExp, args.readOverlap,
                  args.oneSide, args.downsample, args.Hminus, args.Kminus, args.Lminus, args.skipHVR, args.Hr,
               args.Fr)

