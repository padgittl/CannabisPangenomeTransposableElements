import sys, re, os, math, statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import numpy as np
from scipy.stats import pearsonr


########
# MAIN #
########


def readBedFile(bedFile,orientationDict,assemblyLengthDict,assemblyID):
    methylationData = {}
    with open(bedFile,'r') as F:
        for line in F:
            line = line.strip().split()
            scaffoldID = line[0]
            assembly_ID,chrID = line[0].split('.')
            if 'chr' in chrID:
                cytosineStart = line[1]
                cytosineStop = line[2]
                probability = line[10]
                cytosineStart = int(cytosineStart)
                cytosineStop = int(cytosineStop)
                probability = float(probability)
                if assemblyID == assembly_ID:
                    if assemblyID in assemblyLengthDict and chrID in assemblyLengthDict[assemblyID]:
                        chromLength = assemblyLengthDict[assemblyID][chrID]
                    else:
                        print(assemblyID,chrID," not in assemblyLengthDict in readBedFile sub-routine")
                        sys.exit()
                    if assemblyID in orientationDict and chrID in orientationDict[assemblyID]:
                        altBool = orientationDict[assemblyID][chrID]
                        if 'chrY' in chrID:
                            refBool = orientationDict['AH3Mb'][chrID]
                        else:
                            refBool = orientationDict['EH23a'][chrID]
                        if altBool != refBool:
                            newStart = chromLength - cytosineStop
                            newStop = chromLength - cytosineStart
                            cytosineStart = newStart
                            cytosineStop = newStop
                    else:
                        print(assemblyID,chrID," not in orientationDict")
                        sys.exit()
                    if assemblyID not in methylationData:
                        methylationData[assemblyID] = {}
                    if chrID not in methylationData[assemblyID]:
                        methylationData[assemblyID][chrID] = []
                    methylationData[assemblyID][chrID].append((cytosineStart,cytosineStop,probability))
                else:
                    print("assembly IDs don't match")
                    sys.exit()
    return(methylationData)


def calculateMethylationDensity(methylationData,assemblyLengthDict,windowSize,assemblyID,chromID):
    intermediateMethylationDensity = {}
    OUT_METHYLATION = open(assemblyID + '_methylation_densities_window' + str(windowSize) + '_' + chromID + '.txt','w')
    OUT_METHYLATION.write("## window/bin size " + str(windowSize) + " bp\n")
    OUT_METHYLATION.write("## assemblyID\tchromID\tposition\taverageProbability\n")
    for assemblyID in assemblyLengthDict:
        for chromID in assemblyLengthDict[assemblyID]:
            seqLen = assemblyLengthDict[assemblyID][chromID]
            for i in range(0,seqLen,windowSize):
                methylationCount = 0
                position = float(i)
                windowInterval = i + windowSize - 1
                if assemblyID in methylationData and chromID in methylationData[assemblyID]:
                    if assemblyID not in intermediateMethylationDensity:
                        intermediateMethylationDensity[assemblyID] = {}
                    if chromID not in intermediateMethylationDensity[assemblyID]:
                        intermediateMethylationDensity[assemblyID][chromID] = {}
                    if position not in intermediateMethylationDensity[assemblyID][chromID]:
                        intermediateMethylationDensity[assemblyID][chromID][position] = []
                    for cytosineStart,cytosineStop,probability in methylationData[assemblyID][chromID]:
                        if i <= cytosineStart and cytosineStop < windowInterval:
                            methylationCount += 1
                            intermediateMethylationDensity[assemblyID][chromID][position].append(probability)
    methylationDensity = {}
    for assemblyID in intermediateMethylationDensity:
        if assemblyID not in methylationDensity:
            methylationDensity[assemblyID] = {}
        for chromID in intermediateMethylationDensity[assemblyID]:
            if chromID not in methylationDensity[assemblyID]:
                methylationDensity[assemblyID][chromID] = []
            for position in intermediateMethylationDensity[assemblyID][chromID]:
                if len(intermediateMethylationDensity[assemblyID][chromID][position]) > 1:
                    avgProbability = float(sum(intermediateMethylationDensity[assemblyID][chromID][position])) / len(intermediateMethylationDensity[assemblyID][chromID][position])
                else:
                    avgProbability = np.nan
                methylationDensity[assemblyID][chromID].append((position,avgProbability))
                OUT_METHYLATION.write("%s\t%s\t%s\t%s\n" % (assemblyID,chromID,position,avgProbability))
    for assemblyID in methylationDensity:
        for chromID in methylationDensity[assemblyID]:
            methylationDensity[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
    return(methylationDensity)



def readVariantFile(variantFile,assemblyLengthDict,orientationDict,assemblyID):
    variantDict = {}
    with open(variantFile,'r') as F:
        for line in F:
            if "Chrom" not in line:
                chrID,start,stop,assembly_ID = line.strip().split('\t')
                if assembly_ID == assemblyID:
                    start = int(start)
                    stop = int(stop)
                    variantLen = stop - start + 1
                    if assemblyID in assemblyLengthDict and chrID in assemblyLengthDict[assemblyID]:
                        chromLength = assemblyLengthDict[assemblyID][chrID]
                    else:
                        print(assemblyID,chrID," not in assemblyLengthDict in readVariantFile sub-routine")
                        sys.exit()
                    # variants from autosomes/chrX are already oriented relative to EH23a -- only need to be concerned with SVs
                    if assemblyID in orientationDict and chrID in orientationDict[assemblyID]:
                        altBool = orientationDict[assemblyID][chrID]
                        if 'chrY' in chrID:
                            refBool = orientationDict['AH3Mb'][chrID]
                            if altBool != refBool:
                                newStart = chromLength - stop
                                newStop = chromLength - start
                                start = newStart
                                stop = newStop
                    else:
                        print(assemblyID,chrID," not in orientationDict")
                        sys.exit()
                    if assemblyID not in variantDict:
                        variantDict[assemblyID] = {}
                    if chrID not in variantDict[assemblyID]:
                        variantDict[assemblyID][chrID] = []
                    variantDict[assemblyID][chrID].append((start,stop,variantLen))
    return(variantDict)


def removeVariantOverlap(variantDict,assemblyLengthDict):
    overlapFilteredVariantDict = {}
    for assemblyID in variantDict:
        for chromID in variantDict[assemblyID]:
            chromLength = assemblyLengthDict[assemblyID][chromID]
            # print(assemblyID,chromID,chromLength)
            countArray = createCountArray(chromLength)
            # sort on longest sequence
            variantDict[assemblyID][chromID].sort(key=lambda x:x[2], reverse=True)
            if assemblyID not in overlapFilteredVariantDict:
                overlapFilteredVariantDict[assemblyID] = {}
            if chromID not in overlapFilteredVariantDict[assemblyID]:
                overlapFilteredVariantDict[assemblyID][chromID] = []
            for start,stop,variantLen in variantDict[assemblyID][chromID]:
                # print("\t",start,stop,variantLen)
                # based on AH3Mb, chr3, which has a variant that extends beyond end of chrom
                if stop > chromLength:
                    stop = chromLength
                countArray,maxValue = updateCountArray(countArray,start,stop)
                if maxValue == 1:
                    overlapFilteredVariantDict[assemblyID][chromID].append((start,stop))
    for assemblyID in overlapFilteredVariantDict:
        for chromID in overlapFilteredVariantDict[assemblyID]:
            overlapFilteredVariantDict[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
    return(overlapFilteredVariantDict)


def calculateVariantDensity(variantDict,assemblyLengthDict,windowSize,labelKeyword,assemblyID,chromID):
    intermediateVariantDensityDict = {}
    OUT = open(assemblyID + '_' + labelKeyword + '_variant_densities_window' + str(windowSize) + '_' + chromID + '.txt','w')
    OUT.write("assemblyID\tchromID\tposition\tvariantCount\ttotalVariantLength\ttotalVariantPercent\n")
    for assemblyID in assemblyLengthDict:
        for chromID in assemblyLengthDict[assemblyID]:
            seqLen = assemblyLengthDict[assemblyID][chromID]
            # not all chroms are included in the variant file
            if assemblyID in variantDict and chromID in variantDict[assemblyID]:
                if assemblyID not in intermediateVariantDensityDict:
                    intermediateVariantDensityDict[assemblyID] = {}
                if chromID not in intermediateVariantDensityDict[assemblyID]:
                    intermediateVariantDensityDict[assemblyID][chromID] = []
                for i in range(0,seqLen,windowSize):
                    variantCount = 0
                    position = float(i)
                    windowInterval = i + windowSize - 1
                    totalVariantLength = 0
                    if assemblyID in variantDict and chromID in variantDict[assemblyID]:
                        for start,stop in variantDict[assemblyID][chromID]:
                            if i <= start and start < windowInterval:
                                if stop > windowInterval:
                                    stop = windowInterval
                                variantCount += 1
                                variantLen = stop - start + 1
                                totalVariantLength += variantLen
                            if start < i and stop > i and stop < windowInterval:
                                start = i
                                variantCount += 1
                                variantLen = stop - start + 1
                                totalVariantLength += variantLen
                            if start < i and stop > windowInterval:
                                start = i
                                stop = windowInterval
                                variantCount += 1
                                variantLen = stop - start + 1
                                totalVariantLength += variantLen
                    else:
                        print(assemblyID,chromID,"not in variantDict in calculateVariantDensity subroutine")
                        sys.exit()
                    intermediateVariantDensityDict[assemblyID][chromID].append((position,variantCount,totalVariantLength))
    variantDensityDict = {}
    for assemblyID in intermediateVariantDensityDict:
        if assemblyID not in variantDensityDict:
            variantDensityDict[assemblyID] = {}
        for chromID in intermediateVariantDensityDict[assemblyID]:
            if chromID not in variantDensityDict[assemblyID]:
                variantDensityDict[assemblyID][chromID] = []
            seqLen = assemblyLengthDict[assemblyID][chromID]
            lastWindowLength = math.fmod(seqLen, windowSize)
            if lastWindowLength > 0:
                seqLenToLoopThrough = seqLen - lastWindowLength
            else:
                lastWindowLength = windowSize
                seqLenToLoopThrough = seqLen
            intermediateVariantDensityDict[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
            lastPos,lastVarCount,lastTotalVariantLength = intermediateVariantDensityDict[assemblyID][chromID][-1]
            lastTotalVariantPercent = float(lastTotalVariantLength) / lastWindowLength * 100
            variantDensityDict[assemblyID][chromID].append((lastPos,lastVarCount,lastTotalVariantLength,lastTotalVariantPercent))
            for i in range(len(intermediateVariantDensityDict[assemblyID][chromID])-1):
                position,variantCount,totalVariantLength = intermediateVariantDensityDict[assemblyID][chromID][i]
                totalVariantPercent = float(totalVariantLength) / windowSize * 100
                variantDensityDict[assemblyID][chromID].append((position,variantCount,totalVariantLength,totalVariantPercent))
    for assemblyID in variantDensityDict:
        for chromID in variantDensityDict[assemblyID]:
            variantDensityDict[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
            for position,variantCount,totalVariantLength,totalVariantPercent in variantDensityDict[assemblyID][chromID]:
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chromID,position,variantCount,totalVariantLength,totalVariantPercent))
    return(variantDensityDict)


# "if you reverse compliment when "flip" is true, the chromosomes should all be oriented similarly"
# oriented relative to EH23a
# for further analysis, AH3Mb selected as reference orientation for chrY
def readOrientationFile(orientationFile):
    orientationDict = {}
    refAssembly = 'EH23a'
    with open(orientationFile,'r') as F:
        for line in F:
            if 'Sample' not in line:
                firstCol,assemblyID,chromID,flipBool = line.strip().split('\t')
                if assemblyID not in orientationDict:
                    orientationDict[assemblyID] = {}
                if chromID not in orientationDict[assemblyID]:
                    orientationDict[assemblyID][chromID] = flipBool
    return(orientationDict)



def readFastaToCalculateGC(fastaFile,assemblyLengthDict,orientationDict,windowSize,assemblyID,chromID):
    cgData = {}
    lastWindowRemainderDict = {}
    OUT = open(assemblyID + '_gc_info_window' + str(windowSize) + '_' + chromID + '.txt','w')
    OUT.write("## window/bin size " + str(windowSize) + " bp\n")
    OUT.write("## assemblyID\tchromID\tposition\tobserved_vs_expected\tgc_content\n")
    for record in SeqIO.parse(fastaFile,"fasta"):
        assemblyID,chrID = record.id.split('.')
        if 'chr' in chrID:
            record.seq = record.seq.upper()
            total_a_count = record.seq.count('A')
            total_c_count = record.seq.count('C')
            total_g_count = record.seq.count('G')
            total_t_count = record.seq.count('T')
            seqLen = assemblyLengthDict[assemblyID][chrID]
            if assemblyID in orientationDict and chrID in orientationDict[assemblyID]:
                altBool = orientationDict[assemblyID][chrID]
                if 'chrY' in chrID:
                    refBool = orientationDict['AH3Mb'][chrID]
                else:
                    refBool = orientationDict['EH23a'][chrID]
                if altBool != refBool:
                    record.seq = record.seq.reverse_complement()
            else:
                print(assemblyID,chrID," not in orientationDict")
                sys.exit()
            if seqLen != len(record.seq):
                print("seqLen and len(record.seq) are unequal lengths",seqLen,len(record.seq))
                sys.exit()
            # https://academic.oup.com/biostatistics/article/11/3/499/256898
            ###
            if assemblyID not in cgData:
                cgData[assemblyID] = {}
            if chrID not in cgData[assemblyID]:
                cgData[assemblyID][chrID] = []
            ###
            # modulo allows you to get the last length of the sequence that is not divisable by the 1Mb window, like a remainder
            lastWindowLength = math.fmod(seqLen, windowSize)
            if lastWindowLength > 0:
                seqLenToLoopThrough = seqLen - lastWindowLength
            else:
                lastWindowLength = windowSize
                seqLenToLoopThrough = seqLen
            if assemblyID not in lastWindowRemainderDict:
                lastWindowRemainderDict[assemblyID] = {}
            if chrID not in lastWindowRemainderDict[assemblyID]:
                lastWindowRemainderDict[assemblyID][chrID] = lastWindowLength
            ## the variables with the prefix 'final' are for the last window, less than 1Mb
            ## non-inclusive index, which is why 1 is added
            finalSeqWindow = record.seq[int(seqLenToLoopThrough):seqLen+1]
            final_c_count   = finalSeqWindow.count('C')
            final_g_count   = finalSeqWindow.count('G')
            final_cg_count  = finalSeqWindow.count('CG')
            # Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)
            # final_observed_expected = float(final_cg_count * lastWindowLength) / (final_c_count * final_g_count)
            if final_c_count > 0 and final_g_count > 0:
                final_observed_expected = float(final_cg_count * lastWindowLength) / (final_c_count * final_g_count)
            else:
                final_observed_expected = 'nan'
            final_gcContent = gc_fraction(finalSeqWindow)
            cgData[assemblyID][chrID].append((int(seqLenToLoopThrough),final_observed_expected,final_gcContent))
            for i in range(0,int(seqLenToLoopThrough),windowSize):
                windowInterval = i + windowSize - 1
                # biopython is non-inclusive -- need to go 1 up to get last position
                # https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-1/
                seqWindow = record.seq[i:int(windowInterval)+1]
                c_count   = seqWindow.count('C')
                g_count   = seqWindow.count('G')
                cg_count  = seqWindow.count('CG')
                # "DNA methylation enables transposable element-driven genome expansion"
                # https://www.pnas.org/doi/10.1073/pnas.1921719117#sec-3
                # N(CpG)/N, N(CpG)=number of CpG dinucleotides and N=len of seq

                # observed/expected CpG https://www.biostars.org/p/355876/
                # Number of CpG * N / (Number of C * Number of G)
                # https://www.biorxiv.org/content/10.1101/2021.06.01.446659v2.full
                # https://genome.ucsc.edu/ -- UCSC genome browser description -- cites Gardiner-Garden M, Frommer M. CpG islands in vertebrate genomes. J Mol Biol. 1987 Jul 20;196(2):261-82. PMID: 3656447
                if c_count > 0 and g_count > 0:
                    observed_expected = float(cg_count * windowSize) / (c_count * g_count)
                else:
                    observed_expected = 'nan'
                # https://www.bioinformatics.org/sms2/cpg_islands.html
                # "CpG islands are defined as sequence ranges where the Obs/Exp value is greater than 0.6 and the GC content is greater than 50%"
                gcContent = gc_fraction(seqWindow)
                cgData[assemblyID][chrID].append((i,observed_expected,gcContent))
                OUT.write("%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chrID,i,observed_expected,gcContent))
            OUT.write("%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chrID,seqLenToLoopThrough,final_observed_expected,final_gcContent))
    '''
    for assemblyID in cgData:
        for chrID in cgData[assemblyID]:
            cgData[assemblyID][chrID].sort(key=lambda x:x[0], reverse=False)
    '''
    return(cgData,lastWindowRemainderDict)
            

# "You may grep LTR_retrotransposon *.fasta.mod.EDTA.intact.gff3|less -S to get only intact LTRs."
# "For intact LTR information in EDTA results, you can find the identity info from the last column of the gff3 file. You can use this to calculate the age yourself, with T = K / 2µ = (1 - identity) / 2µ, where µ is the mutation rate of your species in the unit of per bp per year. By default, LTR_retriever uses µ = 1.3e-8 from rice."
def readGeneGFF(geneGFF,orientationDict,assemblyID,windowSize,chromID):
    assemblyLengthDict = {}
    geneData = {}
    transcriptCoords = {}
    OUT_chr_lengths = open(assemblyID + '_chromosomeLengths_window' + str(windowSize) + '_' + chromID + '.txt','w')
    OUT_chr_lengths.write("## assemblyID\tchromosomeID\tchromosomeLength\n")
    with open(geneGFF,'r') as F:
        for line in F:
            if '##sequence-region' in line:
                scaffoldInfo = line.strip().split()
                scaffoldID = scaffoldInfo[1]
                scaffoldLength = scaffoldInfo[3]
                # print(scaffoldID,scaffoldLength)
                if 'chr' in scaffoldID:
                    assemblyID,chrID = scaffoldID.split('.')
                    OUT_chr_lengths.write("%s\t%s\t%s\n" % (assemblyID,chrID,scaffoldLength))
                    if assemblyID not in assemblyLengthDict:
                        assemblyLengthDict[assemblyID] = {}
                    if chrID not in assemblyLengthDict[assemblyID]:
                        assemblyLengthDict[assemblyID][chrID] = int(scaffoldLength)
            else:
                if '#' not in line:
                    scaffoldID,source,feature,start,stop,score,strand,frame,attribute = line.strip().split("\t")
                    assemblyID,chrID = scaffoldID.split('.')
                    start = int(start)
                    stop = int(stop)
                    if 'chr' in scaffoldID:
                        if feature == 'mRNA':
                            # print(line)
                            if assemblyID in assemblyLengthDict and chrID in assemblyLengthDict[assemblyID]:
                                chromLength = assemblyLengthDict[assemblyID][chrID]
                            else:
                                print(assemblyID,chrID," not in assemblyLengthDict in readGeneGFF sub-routine")
                                sys.exit()
                            if assemblyID in orientationDict and chrID in orientationDict[assemblyID]:
                                altBool = orientationDict[assemblyID][chrID]
                                if 'chrY' in chrID:
                                    refBool = orientationDict['AH3Mb'][chrID]
                                else:
                                    refBool = orientationDict['EH23a'][chrID]
                                if altBool != refBool:
                                    newStart = chromLength - stop
                                    newStop = chromLength - start
                                    start = newStart
                                    stop = newStop
                            else:
                                print(assemblyID,chrID," not in orientationDict in readGeneGFF sub-routine")
                                sys.exit()
                            getTranscriptID = re.search('ID=(.+);Parent=',attribute)
                            transcriptID = getTranscriptID.group(1)
                            transcriptCoords[transcriptID] = (start,stop)
                            if assemblyID not in geneData:
                                geneData[assemblyID] = {}
                            if chrID not in geneData[assemblyID]:
                                geneData[assemblyID][chrID] = []
                            geneData[assemblyID][chrID].append((transcriptID,strand,start,stop))
    return(assemblyLengthDict,geneData,transcriptCoords)


def calculateGeneDensity(geneData,assemblyLengthDict,windowSize,assemblyID,chromID):
    geneDensity = {}
    OUT_GENES = open(assemblyID + '_gene_densities_window' + str(windowSize) + '_' + chromID + '.txt','w')
    OUT_GENES.write("## window/bin size " + str(windowSize) + " bp\n")
    OUT_GENES.write("## assemblyID\tchromID\tposition\tgeneCount\n")
    for assemblyID in assemblyLengthDict:
        for chromID in assemblyLengthDict[assemblyID]:
            seqLen = assemblyLengthDict[assemblyID][chromID]
            ###
            if assemblyID not in geneDensity:
                geneDensity[assemblyID] = {}
            if chromID not in geneDensity[assemblyID]:
                geneDensity[assemblyID][chromID] = []
            ###
            for i in range(0,seqLen,windowSize):
                geneCount = 0
                position = float(i)
                windowInterval = i + windowSize - 1
                # windowInterval = (i + windowSize - 1) / 2
                if assemblyID in geneData and chromID in geneData[assemblyID]:
                    for transcriptID,strand,transcriptStart,transcriptStop in geneData[assemblyID][chromID]:
                        if i <= transcriptStart and transcriptStart < windowInterval:
                            geneCount += 1
                geneDensity[assemblyID][chromID].append((position,geneCount))
                OUT_GENES.write("%s\t%s\t%s\t%s\n" % (assemblyID,chromID,position,geneCount))
    return(geneDensity)


def readTEAnnoGFF(teAnnoGFF,assemblyLengthDict,orientationDict,assemblyID):
    teFeatureDict = {}
    featureIDs = {}
    ltrFeatureDict = {}
    ltrFeatureIDs = {}
    lambda_value = 6.1e-09
    # ltrDistanceThreshold = 0.9987 # ~100,000 years ago
    # ltrDistanceThreshold = 0.987
    thousandsTimeDivisor = 1000.0
    millionsTimeDivisor  = 1000000.0
    with open(teAnnoGFF,'r') as F:
        for line in F:
            if not line.startswith('#'):
                scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                start = int(start)
                end = int(end)
                teLen = end - start + 1
                assemblyID,chrID = scaffoldID.split('.')
                if 'chr' in chrID:
                    if assemblyID in assemblyLengthDict and chrID in assemblyLengthDict[assemblyID]:
                        chromLength = assemblyLengthDict[assemblyID][chrID]
                    else:
                        print(assemblyID,chrID," not in assemblyLengthDict in readTEAnnoGFF sub-routine")
                        sys.exit()
                    if assemblyID in orientationDict and chrID in orientationDict[assemblyID]:
                        altBool = orientationDict[assemblyID][chrID]
                        if 'chrY' in chrID:
                            refBool = orientationDict['AH3Mb'][chrID]
                        else:
                            refBool = orientationDict['EH23a'][chrID]
                        if altBool != refBool:
                            newStart = chromLength - end
                            newStop = chromLength - start
                            start = newStart
                            end = newStop
                    else:
                        print(assemblyID,chrID," not in orientationDict in readTEAnnoGFF sub-routine")
                        sys.exit()
                
                if 'target_site_duplication' not in feature and 'repeat_region' not in feature and 'long_terminal_repeat' not in feature:
                    # unknown_LTR_retrotransposon
                    if 'LTR_retrotransposon' in feature and 'Gypsy' not in feature and 'Copia' not in feature:
                        feature = feature.replace('LTR_retrotransposon','unknown_LTR_retrotransposon')
                    if 'Gypsy' in feature:
                        feature = feature.replace('Gypsy','Ty3')
                    if 'Copia' in feature:
                        feature = feature.replace('Copia','Ty1')
                    if 'chr' in scaffoldID:
                        if '_retrotransposon' not in feature:
                            # print(attribute)
                            assemblyID,chrID = scaffoldID.split('.')
                            newAttribute = attribute.split(';')
                            suffix,identity = newAttribute[4].split('Identity=')
                            identity = identity.strip()
                            #print(identity)
                            if identity != 'NA':
                                identity = float(identity)
                                dist = 1-identity
                                time = (float(dist) / (2*lambda_value))
                                time_in_thousands = float(time) / thousandsTimeDivisor
                                time_in_millions = float(time) / millionsTimeDivisor
                                if assemblyID not in teFeatureDict:
                                    teFeatureDict[assemblyID] = {}
                                if chrID not in teFeatureDict[assemblyID]:
                                    teFeatureDict[assemblyID][chrID] = {}
                                if feature not in teFeatureDict[assemblyID][chrID]:
                                    teFeatureDict[assemblyID][chrID][feature] = []
                                teFeatureDict[assemblyID][chrID][feature].append((start,end,teLen,time_in_thousands,time_in_millions))
                                if feature not in featureIDs:
                                    featureIDs[feature] = 1
                        
                        else:
                            assemblyID,chrID = scaffoldID.split('.')
                            newAttribute = attribute.split(';')
                            if 'LTRRT_' in newAttribute[0]:
                                suffix,identity = newAttribute[5].split('ltr_identity=')
                                identity = identity.strip()
                                if identity != 'NA':
                                    identity = float(identity)
                                    dist = 1-identity
                                    time = (float(dist) / (2*lambda_value))
                                    time_in_thousands = float(time) / thousandsTimeDivisor
                                    time_in_millions = float(time) / millionsTimeDivisor
                                    if assemblyID not in ltrFeatureDict:
                                        ltrFeatureDict[assemblyID] = {}
                                    if chrID not in ltrFeatureDict[assemblyID]:
                                        ltrFeatureDict[assemblyID][chrID] = {}
                                    if feature not in ltrFeatureDict[assemblyID][chrID]:
                                        ltrFeatureDict[assemblyID][chrID][feature] = []
                                    ltrFeatureDict[assemblyID][chrID][feature].append((start,end,teLen,time_in_thousands,time_in_millions))
                                    if feature not in ltrFeatureIDs:
                                        ltrFeatureIDs[feature] = 1
    return(teFeatureDict,featureIDs,ltrFeatureDict,ltrFeatureIDs)


def createCountArray(seqLen):
    countArray = [0]*seqLen
    return(countArray)


def updateCountArray(countArray,start,stop):
    # the max value corresponds to the maximum position count per coding exon
    maxValue = 0
    for i in range(start,stop):
        countArray[i] += 1
        if maxValue < countArray[i]:
            maxValue = countArray[i]
    return(countArray,maxValue)


def removeOverlap(featureDict,assemblyLengthDict):
    overlapFilteredFeatureDict = {}
    for assemblyID in featureDict:
        for chromID in featureDict[assemblyID]:
            chromLength = assemblyLengthDict[assemblyID][chromID]
            for featureID in featureDict[assemblyID][chromID]:
                countArray = createCountArray(chromLength)
                featureDict[assemblyID][chromID][featureID].sort(key=lambda x:x[2], reverse=True)
                if assemblyID not in overlapFilteredFeatureDict:
                    overlapFilteredFeatureDict[assemblyID] = {}
                if chromID not in overlapFilteredFeatureDict[assemblyID]:
                    overlapFilteredFeatureDict[assemblyID][chromID] = {}
                if featureID not in overlapFilteredFeatureDict[assemblyID][chromID]:
                    overlapFilteredFeatureDict[assemblyID][chromID][featureID] = []
                for teStart,teStop,teLen,time_in_thousands,time_in_millions in featureDict[assemblyID][chromID][featureID]:
                    countArray,maxValue = updateCountArray(countArray,teStart,teStop)
                    if maxValue == 1:
                        overlapFilteredFeatureDict[assemblyID][chromID][featureID].append((teStart,teStop,teLen,time_in_thousands,time_in_millions))
    for assemblyID in overlapFilteredFeatureDict:
        for chromID in overlapFilteredFeatureDict[assemblyID]:
            for featureID in overlapFilteredFeatureDict[assemblyID][chromID]:
                overlapFilteredFeatureDict[assemblyID][chromID][featureID].sort(key=lambda x:x[0], reverse=False)
    return(overlapFilteredFeatureDict)

'''
CACTA_TIR_transposon
Copia_LTR_retrotransposon
Gypsy_LTR_retrotransposon
LTR_retrotransposon
Mutator_TIR_transposon
PIF_Harbinger_TIR_transposon
Tc1_Mariner_TIR_transposon
hAT_TIR_transposon
helitron
'''

def createColorMapDict(featureIDs,featureType):
    colorMap = {}
    featureIDList = list(featureIDs.keys())

    teColors = {'CACTA_TIR_transposon':'#000000','Mutator_TIR_transposon':'#E69F00','PIF_Harbinger_TIR_transposon':'#0072B2','Tc1_Mariner_TIR_transposon':'#009E73','hAT_TIR_transposon':'#F0E442','helitron':'#CC79A7'}

    ltrColors = {'Ty1_LTR_retrotransposon':'#542788', 'unknown_LTR_retrotransposon':'#998ec3', 'Ty3_LTR_retrotransposon':'#b35806'}
    if featureType == 'TEs':
        for featureID in featureIDList:
            if featureID not in colorMap:
                colorMap[featureID] = teColors[featureID]
    if featureType == 'LTRs':
        for featureID in featureIDList:
            if featureID not in colorMap:
                colorMap[featureID] = ltrColors[featureID]
    return(colorMap)


# windows -->
# 0-5, 6-10, 11-15, 16-20, 21-25
def calculateTEDensity(assemblyLengthDict,teFeatureDict,assemblyID,featureType,windowSize,chromID):
    intermediateCompiledTimes = {}
    intermediateDensityData = {}
    OUT_TIMES =	open(assemblyID + '_' + featureType + '_times_window' + str(windowSize) + '_' + chromID + '.txt','w')
    OUT_TIMES.write("## window/bin size " + str(windowSize) + " bp\n")
    OUT_TIMES.write("## assemblyID\tchromID\tfeatureID\ttimeID\tchromStartPosition\tspecificTECount\tteCount\ttotalTELenPerWindow\tpercentWindow\taverageTime\tspecificTimeTEPercent\n")
    OUT_TIMES.write("## Time units are in millions of years\n")
    OUT_DENSITIES = open(assemblyID + '_' + featureType + '_densities_window' + str(windowSize) + '_' + chromID + '.txt','w')
    OUT_DENSITIES.write("## window/bin size " + str(windowSize) + " bp\n")
    OUT_DENSITIES.write("## assemblyID\tchromID\tfeatureID\tchromStartPosition\tteCount\ttotalTELenPerWindow\tpercentWindow\n")
    for assemblyID in assemblyLengthDict:
        for chromID in assemblyLengthDict[assemblyID]:
            seqLen = assemblyLengthDict[assemblyID][chromID]
            if assemblyID in teFeatureDict and chromID in teFeatureDict[assemblyID]:
                if assemblyID not in intermediateDensityData:
                    intermediateDensityData[assemblyID] = {}
                if chromID not in intermediateDensityData[assemblyID]:
                    intermediateDensityData[assemblyID][chromID] = {}
                if assemblyID not in intermediateCompiledTimes:
                    intermediateCompiledTimes[assemblyID] = {}
                if chromID not in intermediateCompiledTimes[assemblyID]:
                    intermediateCompiledTimes[assemblyID][chromID] = {}
                for featureID in teFeatureDict[assemblyID][chromID]:
                    if featureID not in intermediateDensityData[assemblyID][chromID]:
                        intermediateDensityData[assemblyID][chromID][featureID] = []
                    if featureID not in intermediateCompiledTimes[assemblyID][chromID]:
                        intermediateCompiledTimes[assemblyID][chromID][featureID] = {}
                    # at this point, getting counts for different TEs, not calculating percent, so the last window length doesn't matter yet
                    # looping through all windows for every TE feature
                    for i in range(0,seqLen,windowSize):
                        totalTELen = 0
                        teCount = 0
                        position = float(i)
                        windowInterval = i + windowSize - 1
                        
                        teCount0_5 = 0
                        teLen0_5 = 0
                        teList0_5 = []
                        
                        teCount6_10 = 0
                        teLen6_10 = 0
                        teList6_10 = []
                        
                        teCount11_15 = 0
                        teLen11_15 = 0
                        teList11_15 = []
                        
                        teCount16_20 = 0
                        teLen16_20 = 0
                        teList16_20 = []
                        
                        teCount21_25 = 0
                        teLen21_25 = 0
                        teList21_25 = []
                        
                        for teStart,teStop,teLen,time_in_thousands,time_in_millions in teFeatureDict[assemblyID][chromID][featureID]:
                            # require that full feature is contained within window, otherwise this causes a problem with overlapping sequences
                            # in this example, the TE starts inside the window and either terminates inside the window or outside of it
                            if i <= teStart and teStart < windowInterval:
                                if teStop > windowInterval:
                                    teStop = windowInterval
                                teCount += 1
                                newTELen = teStop - teStart + 1
                                totalTELen += newTELen
                                # print(i,windowInterval,time_in_thousands,time_in_millions)
                                if time_in_millions < 6:
                                    teCount0_5 += 1
                                    teLen0_5 += newTELen
                                    teList0_5.append(time_in_millions)
                                if time_in_millions >= 6 and time_in_millions < 11:
                                    teCount6_10 += 1
                                    teLen6_10 += newTELen
                                    teList6_10.append(time_in_millions)
                                if time_in_millions >= 11 and time_in_millions < 16:
                                    teCount11_15 += 1
                                    teLen11_15 += newTELen
                                    teList11_15.append(time_in_millions)
                                if time_in_millions >= 16 and time_in_millions < 21:
                                    teCount16_20 += 1
                                    teLen16_20 += newTELen
                                    teList16_20.append(time_in_millions)
                                if time_in_millions >= 21 and time_in_millions < 26:
                                    teCount21_25 += 1
                                    teLen21_25 += newTELen
                                    teList21_25.append(time_in_millions)

                            # in this example, teStart comes before the window starts, but terminates inside the window
                            if teStart < i and teStop > i and teStop < windowInterval:
                                teStart = i
                                teCount += 1
                                newTELen = teStop - teStart
                                totalTELen += newTELen
                                if time_in_millions < 6:
                                    teCount0_5 += 1
                                    teLen0_5 += newTELen
                                    teList0_5.append(time_in_millions)
                                if time_in_millions >= 6 and time_in_millions < 11:
                                    teCount6_10 += 1
                                    teLen6_10 += newTELen
                                    teList6_10.append(time_in_millions)
                                if time_in_millions >= 11 and time_in_millions < 16:
                                    teCount11_15 += 1
                                    teLen11_15 += newTELen
                                    teList11_15.append(time_in_millions)
                                if time_in_millions >= 16 and time_in_millions < 21:
                                    teCount16_20 += 1
                                    teLen16_20 += newTELen
                                    teList16_20.append(time_in_millions)
                                if time_in_millions >= 21 and time_in_millions < 26:
                                    teCount21_25 += 1
                                    teLen21_25 += newTELen
                                    teList21_25.append(time_in_millions)
                                    
                            # this example covers the whole window
                            if teStart < i and teStop > windowInterval:
                                teStart = i
                                teStop = windowInterval
                                newTELen = teStop - teStart + 1
                                totalTELen += newTELen
                                teCount += 1
                                if time_in_millions < 6:
                                    teCount0_5 += 1
                                    teLen0_5 += newTELen
                                    teList0_5.append(time_in_millions)
                                if time_in_millions >= 6 and time_in_millions < 11:
                                    teCount6_10 += 1
                                    teLen6_10 += newTELen
                                    teList6_10.append(time_in_millions)
                                if time_in_millions >= 11 and time_in_millions < 16:
                                    teCount11_15 += 1
                                    teLen11_15 += newTELen
                                    teList11_15.append(time_in_millions)
                                if time_in_millions >= 16 and time_in_millions < 21:
                                    teCount16_20 += 1
                                    teLen16_20 += newTELen
                                    teList16_20.append(time_in_millions)
                                if time_in_millions >= 21 and time_in_millions < 26:
                                    teCount21_25 += 1
                                    teLen21_25 += newTELen
                                    teList21_25.append(time_in_millions)
                        # at this point, the time window for a TE has been counted in the above for-loop
                        if len(teList0_5) > 0:
                            avg0_5 = float(sum(teList0_5)) / len(teList0_5)
                        else:
                            avg0_5 = 0
                        if len(teList6_10) > 0:
                            avg6_10 = float(sum(teList6_10)) / len(teList6_10)
                        else:
                            avg6_10 = 0
                        if len(teList11_15) > 0:
                            avg11_15 = float(sum(teList11_15)) / len(teList11_15)
                        else:
                            avg11_15 = 0
                        if len(teList16_20) > 0:
                            avg16_20 = float(sum(teList16_20)) / len(teList16_20)
                        else:
                            avg16_20 = 0
                        if len(teList21_25) > 0:
                            avg21_25 = float(sum(teList21_25)) / len(teList21_25)
                        else:
                            avg21_25 = 0
                        # print(assemblyID,chromID,featureID,position,teCount)
                        intermediateDensityData[assemblyID][chromID][featureID].append((position,teCount,totalTELen))
                        if '0_5' not in intermediateCompiledTimes[assemblyID][chromID][featureID]:
                            intermediateCompiledTimes[assemblyID][chromID][featureID]['0_5'] = []
                        intermediateCompiledTimes[assemblyID][chromID][featureID]['0_5'].append((position,teCount0_5,teCount,teLen0_5,totalTELen,avg0_5))
                        if '6_10' not in intermediateCompiledTimes[assemblyID][chromID][featureID]:
                            intermediateCompiledTimes[assemblyID][chromID][featureID]['6_10'] = []
                        intermediateCompiledTimes[assemblyID][chromID][featureID]['6_10'].append((position,teCount6_10,teCount,teLen6_10,totalTELen,avg6_10))
                        if '11_15' not in intermediateCompiledTimes[assemblyID][chromID][featureID]:
                            intermediateCompiledTimes[assemblyID][chromID][featureID]['11_15'] = []
                        intermediateCompiledTimes[assemblyID][chromID][featureID]['11_15'].append((position,teCount11_15,teCount,teLen11_15,totalTELen,avg11_15))
                        if '16_20' not in intermediateCompiledTimes[assemblyID][chromID][featureID]:
                            intermediateCompiledTimes[assemblyID][chromID][featureID]['16_20'] = []
                        intermediateCompiledTimes[assemblyID][chromID][featureID]['16_20'].append((position,teCount16_20,teCount,teLen16_20,totalTELen,avg16_20))
                        if '21_25' not in intermediateCompiledTimes[assemblyID][chromID][featureID]:
                            intermediateCompiledTimes[assemblyID][chromID][featureID]['21_25'] = []
                        intermediateCompiledTimes[assemblyID][chromID][featureID]['21_25'].append((position,teCount21_25,teCount,teLen21_25,totalTELen,avg21_25))
    compiledTimes = {}
    densityData = {}
    for assemblyID in intermediateDensityData:
        ###
        if assemblyID not in densityData:
            densityData[assemblyID] = {}
        ###
        for chromID in intermediateDensityData[assemblyID]:
            seqLen = assemblyLengthDict[assemblyID][chromID]
            lastWindowLength = math.fmod(seqLen, windowSize)
            if lastWindowLength > 0:
                seqLenToLoopThrough = seqLen - lastWindowLength
            else:
                lastWindowLength = windowSize
                seqLenToLoopThrough = seqLen
            ###
            if chromID not in densityData[assemblyID]:
                densityData[assemblyID][chromID] = {}
            ###
            for featureID in intermediateDensityData[assemblyID][chromID]:
                ###
                if featureID not in densityData[assemblyID][chromID]:
                    densityData[assemblyID][chromID][featureID] = []
                ###
                intermediateDensityData[assemblyID][chromID][featureID].sort(key=lambda x:x[0], reverse=False)
                lastPos,lastTECount,lastTotalTELen = intermediateDensityData[assemblyID][chromID][featureID][-1]
                lastTotalTEPercent = float(lastTotalTELen) / lastWindowLength * 100
                densityData[assemblyID][chromID][featureID].append((lastPos,lastTECount,lastTotalTELen,lastTotalTEPercent))
                for i in range(len(intermediateDensityData[assemblyID][chromID][featureID])-1):
                    position,teCount,totalTELen = intermediateDensityData[assemblyID][chromID][featureID][i]
                    totalTEPercent = float(totalTELen) / windowSize * 100
                    densityData[assemblyID][chromID][featureID].append((position,teCount,totalTELen,totalTEPercent))
                    OUT_DENSITIES.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chromID,featureID,position,teCount,totalTELen,totalTEPercent))
                OUT_DENSITIES.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chromID,featureID,lastPos,lastTECount,lastTotalTELen,lastTotalTEPercent))
    for assemblyID in densityData:
        for chromID in densityData[assemblyID]:
            for featureID in densityData[assemblyID][chromID]:
                densityData[assemblyID][chromID][featureID].sort(key=lambda x:x[0], reverse=False)
    for assemblyID in intermediateCompiledTimes:
        ###
        if assemblyID not in compiledTimes:
            compiledTimes[assemblyID] = {}
        ###
        for chromID in intermediateCompiledTimes[assemblyID]:
            seqLen = assemblyLengthDict[assemblyID][chromID]
            lastWindowLength = math.fmod(seqLen, windowSize)
            if lastWindowLength > 0:
                seqLenToLoopThrough = seqLen - lastWindowLength
            else:
                lastWindowLength = windowSize
                seqLenToLoopThrough = seqLen
            ###
            if chromID not in compiledTimes[assemblyID]:
                compiledTimes[assemblyID][chromID] = {}
            ###
            for featureID in intermediateCompiledTimes[assemblyID][chromID]:
                ###
                if featureID not in compiledTimes[assemblyID][chromID]:
                    compiledTimes[assemblyID][chromID][featureID] = {}
                ###
                for timeID in intermediateCompiledTimes[assemblyID][chromID][featureID]:
                    ###
                    if timeID not in compiledTimes[assemblyID][chromID][featureID]:
                        compiledTimes[assemblyID][chromID][featureID][timeID] = []
                    ###
                    # intermediateCompiledTimes[assemblyID][chromID][featureID]['21_25'].append((position,teCount,teLen21_25,totalTELen,avg21_25))
                    intermediateCompiledTimes[assemblyID][chromID][featureID][timeID].sort(key=lambda x:x[0], reverse=False)
                    lastPos,lastSpecificTECount,lastTECount,lastTELenForSpecificTime,lastTotalTELen,lastAvgTime = intermediateCompiledTimes[assemblyID][chromID][featureID][timeID][-1]
                    # the TE percent is relative to the length of the window, not the length of all of the differently-aged TEs combined
                    lastTotalTEPercent = float(lastTotalTELen) / lastWindowLength * 100
                    lastSpecificTimeTEPercent = float(lastTELenForSpecificTime) / lastWindowLength * 100
                    # print(assemblyID,chromID,featureID,timeID,lastTotalTEPercent,lastSpecificTimeTEPercent)
                    for i in range(len(intermediateCompiledTimes[assemblyID][chromID][featureID][timeID])-1):
                        position,specificTECount,teCount,teLenForSpecificTime,totalTELen,avgTime = intermediateCompiledTimes[assemblyID][chromID][featureID][timeID][i]
                        totalTEPercent = float(totalTELen) / windowSize * 100
                        specificTimeTEPercent = float(teLenForSpecificTime) / windowSize * 100
                        # print(assemblyID,chromID,featureID,timeID,position,teCount,totalTELen,totalTEPercent,avgTime,specificTimeTEPercent)
                        compiledTimes[assemblyID][chromID][featureID][timeID].append((position,specificTECount,teCount,totalTELen,totalTEPercent,avgTime,specificTimeTEPercent))
                        OUT_TIMES.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chromID,featureID,timeID,position,specificTECount,teCount,totalTELen,totalTEPercent,avgTime,specificTimeTEPercent))
                    compiledTimes[assemblyID][chromID][featureID][timeID].append((lastPos,lastSpecificTECount,lastTECount,lastTotalTELen,lastTotalTEPercent,lastAvgTime,lastSpecificTimeTEPercent))
                    OUT_TIMES.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chromID,featureID,timeID,lastPos,lastSpecificTECount,lastTECount,lastTotalTELen,lastTotalTEPercent,lastAvgTime,lastSpecificTimeTEPercent))
                    # NLv1b   chr7    Mutator_TIR_transposon  0_5     67090000.0      2       991     99.1    3.729508196721315       48856.5
                    # NLv1b   chr7    Mutator_TIR_transposon  0_5     66846000.0      2       991     99.1    3.729508196721315       38368.3
                    # NLv1b   chr7    Mutator_TIR_transposon  0_5     371000.0        2       899     89.9    2.9918032786885274      30888.4
                    # NLv1b   chr7    Mutator_TIR_transposon  0_5     67566000.0      2       899     89.9    2.9918032786885274      28232.100000000002
    '''
    for assemblyID in compiledTimes:
        for chromID in compiledTimes[assemblyID]:
            for featureID in compiledTimes[assemblyID][chromID]:
                for timeID in compiledTimes[assemblyID][chromID][featureID]:
                    compiledTimes[assemblyID][chromID][featureID][timeID].sort(key=lambda x:x[0], reverse=False)
    '''
    return(densityData,compiledTimes)
                    


usage = "Usage: " + sys.argv[0] + " <gene gff> <TE annotation gff> <methylation bed file> <fasta file> <orientation file> <dup variant file> <inversion variant file> <translocation variant file> <assembly ID> <chromosome ID> <window size, e.g. 1000000>\n"
if len(sys.argv) != 12:
    print(usage)
    sys.exit()

geneGFF = sys.argv[1]
teAnnoGFF = sys.argv[2]
bedFile = sys.argv[3]
fastaFile = sys.argv[4]
orientationFile = sys.argv[5]
dupVariantFile = sys.argv[6]
invVariantFile = sys.argv[7]
transVariantFile = sys.argv[8]
assemblyID = sys.argv[9]
chromID = sys.argv[10]
windowSize = sys.argv[11]

# windowSize = 1000000
windowSize = int(windowSize)
obs_exp_threshold = 0.6

orientationDict = readOrientationFile(orientationFile)
assemblyLengthDict,geneData,transcriptCoords = readGeneGFF(geneGFF,orientationDict,assemblyID,windowSize,chromID)

methylationData = readBedFile(bedFile,orientationDict,assemblyLengthDict,assemblyID)
methylationDensity = calculateMethylationDensity(methylationData,assemblyLengthDict,windowSize,assemblyID,chromID)

dupVariantDict = readVariantFile(dupVariantFile,assemblyLengthDict,orientationDict,assemblyID)
invVariantDict = readVariantFile(invVariantFile,assemblyLengthDict,orientationDict,assemblyID)
transVariantDict = readVariantFile(transVariantFile,assemblyLengthDict,orientationDict,assemblyID)

dup_overlapFilteredVariantDict = removeVariantOverlap(dupVariantDict,assemblyLengthDict)
# print here
dup_variantDensityDict = calculateVariantDensity(dup_overlapFilteredVariantDict,assemblyLengthDict,windowSize,'dups',assemblyID,chromID)

inv_overlapFilteredVariantDict = removeVariantOverlap(invVariantDict,assemblyLengthDict)
# print here
inv_variantDensityDict = calculateVariantDensity(inv_overlapFilteredVariantDict,assemblyLengthDict,windowSize,'inv',assemblyID,chromID)

trans_overlapFilteredVariantDict = removeVariantOverlap(transVariantDict,assemblyLengthDict)
# print here
trans_variantDensityDict = calculateVariantDensity(trans_overlapFilteredVariantDict,assemblyLengthDict,windowSize,'transloc',assemblyID,chromID)

# print here
geneDensity = calculateGeneDensity(geneData,assemblyLengthDict,windowSize,assemblyID,chromID)

# print here
cgData,lastWindowRemainderDict = readFastaToCalculateGC(fastaFile,assemblyLengthDict,orientationDict,windowSize,assemblyID,chromID)

teFeatureDict,featureIDs,ltrFeatureDict,ltrFeatureIDs = readTEAnnoGFF(teAnnoGFF,assemblyLengthDict,orientationDict,assemblyID)

overlapFilteredTEFeatureDict = removeOverlap(teFeatureDict,assemblyLengthDict)
overlapFilteredLTRFeatureDict = removeOverlap(ltrFeatureDict,assemblyLengthDict)

# print here
teDensityData,compiledTimes = calculateTEDensity(assemblyLengthDict,overlapFilteredTEFeatureDict,assemblyID,'TEs',windowSize,chromID)
ltrDensityData,compiledLTRTimes = calculateTEDensity(assemblyLengthDict,overlapFilteredLTRFeatureDict,assemblyID,'LTRs',windowSize,chromID)





