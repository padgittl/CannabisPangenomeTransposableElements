import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

########
# MAIN #
########

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6265445/
# "A sequence region is termed solo LTR if (i) it is annotated by an LTR region without any internal regions located within 300 bp adjacent to the target region; (ii) no nearby (the adjacent four annotation entries) sequence regions were annotated by the same LTR-RT entry and (iii) the length of the alignment hit accounts for at least 80% of the length of the solo LTR candidate."

# 1. min 300 bp distance between candidate soloLTR and intLTR
# 2. adjacent four annotation entries should not be from the same LTR
# 3. similarity greater than 0.8 perc iden


def readCannabisGeneData(gffFile):
    chrLengthsDict = {}
    with open(gffFile,'r') as F:
        for line in F:
            if '##sequence-region' in line:
                # ##sequence-region AH3Ma.100 1 125601
                prefix,scaffoldID,chromStart,chromosomeLength = line.strip().split()
                if 'chr' in scaffoldID:
                    assemblyID,chrID = scaffoldID.split('.')
                    chromosomeLength = int(chromosomeLength)
                    # print(assemblyID,chrID,chromosomeLength)
                    if assemblyID not in chrLengthsDict:
                        chrLengthsDict[assemblyID] = {}
                    if chrID not in chrLengthsDict[assemblyID]:
                        chrLengthsDict[assemblyID][chrID] = chromosomeLength
    return(chrLengthsDict)


def readMaizeGeneData(gffFile):
    chrLengthsDict = {}
    basename = os.path.basename(gffFile)
    assemblyID,suffixID = basename.strip().split('-REFERENCE')
    if assemblyID not in chrLengthsDict:
        chrLengthsDict[assemblyID] = {}
    with open(gffFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                if feature == 'chromosome':
                    chromosomeLength = int(end)
                    if scaffoldID not in chrLengthsDict[assemblyID]:
                        chrLengthsDict[assemblyID][scaffoldID] = chromosomeLength
    return(chrLengthsDict)


# https://www.nature.com/articles/s41467-017-02546-5#Sec7
def readTEAnnoGFF(gffFile,specificStrand,specificAssemblyID,organismID):
    minAlignmentScore = 300
    minAlignmentLength = 100
    minPercentIdentity = 0.8
    candidateSoloLTRCoordDict = {}
    internalFeatureCoordDict = {}
    internalSeqLengthDict = {}
    intactLTRCoordDict = {}
    countDict = {}
    OUT = open(specificAssemblyID + '_unfilteredSoloLTRs.tsv','w')
    OUT.write("assemblyID\tchrID\tltrType\tteID\tuniqueTEID\tstart\tend\talignLen\tstrand\tscore\tpercIden\n")
    with open(gffFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                if organismID == 'cannabis':
                    assemblyID,chrID = scaffoldID.split('.')
                elif organismID == 'maize':
                    basename = os.path.basename(gffFile)
                    assemblyID,suffixID = basename.strip().split('-REFERENCE')
                    chrID = scaffoldID
                start = int(start)
                end = int(end)
                alignLen = end - start + 1
                if 'chr' in scaffoldID:
                    # this is counting ltr features and internal sequence features
                    if 'Gypsy' not in countDict:
                        countDict['Gypsy'] = {}
                    if assemblyID not in countDict['Gypsy']:
                        countDict['Gypsy'][assemblyID] = {}
                    if chrID not in countDict['Gypsy'][assemblyID]:
                        countDict['Gypsy'][assemblyID][chrID] = {}
                    if 'INT' not in countDict['Gypsy'][assemblyID][chrID]:
                        countDict['Gypsy'][assemblyID][chrID]['INT'] = 0
                    if 'LTR' not in countDict['Gypsy'][assemblyID][chrID]:
                        countDict['Gypsy'][assemblyID][chrID]['LTR'] = 0
                    ### 
                    ### 
                    if 'Copia' not in countDict:
                        countDict['Copia'] = {}
                    if assemblyID not in countDict['Copia']:
                        countDict['Copia'][assemblyID] = {}
                    if chrID not in countDict['Copia'][assemblyID]:
                        countDict['Copia'][assemblyID][chrID] = {}
                    if 'INT' not in countDict['Copia'][assemblyID][chrID]:
                        countDict['Copia'][assemblyID][chrID]['INT'] = 0
                    if 'LTR' not in countDict['Copia'][assemblyID][chrID]:
                        countDict['Copia'][assemblyID][chrID]['LTR'] = 0
                            
                    if 'LTR_retrotransposon' in feature:
                        newAttr = attribute.split(';')
                        method = newAttr[-1]
                        # Method=structural are intact -- not looking for these examples here
                        if method == 'Method=homology':
                            # print(attribute,strand,specificStrand)
                            if strand == specificStrand:
                                # print(strand,specificStrand)
                                prefix,percIden = newAttr[-2].split('Identity=')
                                percIden = float(percIden)
                                score = int(score)
                                teID = newAttr[1]
                                prefix,teID = teID.split('Name=')
                                prefix,uniqueTEID = newAttr[0].split('ID=')
                                # print(teID,uniqueTEID)
                                # below are the thresholds specifically for solo-LTRs
                                if score > minAlignmentScore:
                                    if alignLen >= minAlignmentLength:
                                        if percIden >= minPercentIdentity:
                                            # not an internal sequence, but may not explicitly include 'LTR' in name
                                            if '_INT' not in teID:
                                                if '_LTR' not in teID:
                                                    teID = teID + '_LTR'
                                                    # print(teID)
                                                    # if score > minAlignmentScore:
                                                    # $end-$start<=100
                                                    # "at least 100 bp in length"
                                                    # if alignLen >= minAlignmentLength:
                                                    # https://github.com/oushujun/LTR_retriever/blob/master/bin/solo_finder.pl
                                                    # "$keep=0 if $loc2=~/int/i; #skip this line if encounter internal regions"
                                                    # 80% of the library LTR  
                                                    # if percIden >= minPercentIdentity:
                                                if 'Gypsy' in feature:
                                                    OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chrID,'Ty3',teID,uniqueTEID,start,end,alignLen,strand,score,percIden))
                                                    countDict['Gypsy'][assemblyID][chrID]['LTR'] += 1

                                                    if 'Gypsy' not in candidateSoloLTRCoordDict:
                                                        candidateSoloLTRCoordDict['Gypsy'] = {}
                                                    if chrID not in candidateSoloLTRCoordDict['Gypsy']:
                                                        candidateSoloLTRCoordDict['Gypsy'][chrID] = {}
                                                    if assemblyID not in candidateSoloLTRCoordDict['Gypsy'][chrID]:
                                                        candidateSoloLTRCoordDict['Gypsy'][chrID][assemblyID] = []
                                                    candidateSoloLTRCoordDict['Gypsy'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID,strand))
                                                # print(start,end,alignLen,teID,uniqueTEID)
                                                elif 'Copia' in feature:
                                                    OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chrID,'Ty1',teID,uniqueTEID,start,end,alignLen,strand,score,percIden))
                                                    countDict['Copia'][assemblyID][chrID]['LTR'] += 1
                                                    #print(assemblyID,chrID,countDict['Copia'][assemblyID][chrID]['LTR'])
                                                    if 'Copia' not in candidateSoloLTRCoordDict:
                                                        candidateSoloLTRCoordDict['Copia'] = {}
                                                    if chrID not in candidateSoloLTRCoordDict['Copia']:
                                                        candidateSoloLTRCoordDict['Copia'][chrID] = {}
                                                    if assemblyID not in candidateSoloLTRCoordDict['Copia'][chrID]:
                                                        candidateSoloLTRCoordDict['Copia'][chrID][assemblyID] = []
                                                    candidateSoloLTRCoordDict['Copia'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID,strand))
                                                    '''
                                                    # check above where if-statement applies --> to Method=homology
                                                    # https://github.com/oushujun/LTR_retriever/blob/master/bin/solo_finder.pl
                                                    # "if there is int exists in this working window, require it has at least 300bp distance from a solo-LTR"
                                                    '''
                                            if '_INT' in teID:
                                                # print(teID)
                                                if 'Gypsy' in feature:
                                                    countDict['Gypsy'][assemblyID][chrID]['INT'] += 1
                                                    # print(assemblyID,chrID,countDict['Gypsy'][assemblyID][chrID]['INT'])
                                                    if 'Gypsy' not in internalFeatureCoordDict:
                                                        internalFeatureCoordDict['Gypsy'] = {}
                                                    if chrID not in internalFeatureCoordDict['Gypsy']:
                                                        internalFeatureCoordDict['Gypsy'][chrID] = {}
                                                    if assemblyID not in internalFeatureCoordDict['Gypsy'][chrID]:
                                                        internalFeatureCoordDict['Gypsy'][chrID][assemblyID] = []
                                                    internalFeatureCoordDict['Gypsy'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID,strand))
                                                    # print(chrID,assemblyID,teID,uniqueTEID)
                                                    if 'Gypsy' not in internalSeqLengthDict:
                                                        internalSeqLengthDict['Gypsy'] = []
                                                    internalSeqLengthDict['Gypsy'].append(alignLen)
                                                elif 'Copia' in feature:
                                                    countDict['Copia'][assemblyID][chrID]['INT'] += 1
                                                    # print(assemblyID,chrID,countDict['Copia'][assemblyID][chrID]['INT'])
                                                    if 'Copia' not in internalFeatureCoordDict:
                                                        internalFeatureCoordDict['Copia'] = {}
                                                    if chrID not in internalFeatureCoordDict['Copia']:
                                                        internalFeatureCoordDict['Copia'][chrID] = {}
                                                    if assemblyID not in internalFeatureCoordDict['Copia'][chrID]:
                                                        internalFeatureCoordDict['Copia'][chrID][assemblyID] = []
                                                    internalFeatureCoordDict['Copia'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID,strand))
                                                    if 'Copia' not in internalSeqLengthDict:
                                                        internalSeqLengthDict['Copia'] = []
                                                    internalSeqLengthDict['Copia'].append(alignLen)
                        else:
                            method = newAttr[-3]
                            if method == 'Method=structural':
                                if strand == specificStrand:
                                    # ID=LTRRT_1177;Parent=repeat_region_1177;Name=TE_00004144;Classification=LTR/Copia;Sequence_ontology=SO:0002264;ltr_identity=1.0000;Method=structural;motif=TACA;tsd=AATAT
                                    prefix,percIden = newAttr[-4].split('ltr_identity=')
                                    percIden = float(percIden)
                                    # method=structural, the score is given as '.'
                                    teID = newAttr[2]
                                    prefix,teID = teID.split('Name=')
                                    prefix,uniqueTEID = newAttr[0].split('ID=')
                                    # print(teID,uniqueTEID)
                                    if alignLen >= minAlignmentLength:
                                        if percIden >= minPercentIdentity:
                                            if 'Gypsy' in feature:
                                                if 'Gypsy' not in intactLTRCoordDict:
                                                    intactLTRCoordDict['Gypsy'] = {}
                                                if chrID not in intactLTRCoordDict['Gypsy']:
                                                    intactLTRCoordDict['Gypsy'][chrID] = {}
                                                if assemblyID not in intactLTRCoordDict['Gypsy'][chrID]:
                                                    intactLTRCoordDict['Gypsy'][chrID][assemblyID] = []
                                                intactLTRCoordDict['Gypsy'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID,strand))
                                            elif 'Copia' in feature:
                                                if 'Copia' not in intactLTRCoordDict:
                                                    intactLTRCoordDict['Copia'] = {}
                                                if chrID not in intactLTRCoordDict['Copia']:
                                                    intactLTRCoordDict['Copia'][chrID] = {}
                                                if assemblyID not in intactLTRCoordDict['Copia'][chrID]:
                                                    intactLTRCoordDict['Copia'][chrID][assemblyID] = []
                                                intactLTRCoordDict['Copia'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID,strand))
                                                
    for ltrType in countDict:
        for assemblyID in countDict[ltrType]:
            for chrID in countDict[ltrType][assemblyID]:
                intCount = countDict[ltrType][assemblyID][chrID]['INT']
                ltrCount = countDict[ltrType][assemblyID][chrID]['LTR']
                # print(ltrType,assemblyID,chrID,intCount,ltrCount)
    return(candidateSoloLTRCoordDict,internalFeatureCoordDict,internalSeqLengthDict,countDict,intactLTRCoordDict)


# candidateSoloLTRCoordDict['Gypsy'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID))
# # 2. adjacent four annotation entries should not be from the same LTR
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6265445/
# "no nearby (the adjacent four annotation entries) sequence regions were annotated by the same LTR-RT entry"
# double-checked this subroutine on 02/03/2023
def removeMatchingAdjacentEntries(candidateSoloLTRCoordDict):
    # these are examples to keep
    # 5' dir
    upstreamMatchDict = {}
    # 3' dir
    downstreamMatchDict = {}
    for ltrType in candidateSoloLTRCoordDict:
        # don't necessarily know enough about 'Other' LTR type to do this comparison
        if ltrType != 'Other':
            for chrID in candidateSoloLTRCoordDict[ltrType]:
                for assemblyID in candidateSoloLTRCoordDict[ltrType][chrID]:
                    # sort by start position, ascending
                    candidateSoloLTRCoordDict[ltrType][chrID][assemblyID].sort(key = lambda x:x[0], reverse=False)
                    # the "-2" is to allow i+1 and i+2 searches
                    # for i in range(len(candidateSoloLTRCoordDict[ltrType][chrID][assemblyID])-2):
                    for i in range(len(candidateSoloLTRCoordDict[ltrType][chrID][assemblyID])-4):
                        # 4 upstream indexes
                        start_n1,end_n1,alignLen_n1,teID_n1,uniqueTEID_n1,strand_n1 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i-1]
                        start_n2,end_n2,alignLen_n2,teID_n2,uniqueTEID_n2,strand_n2 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i-2]
                        start_n3,end_n3,alignLen_n3,teID_n3,uniqueTEID_n3,strand_n3 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i-3]
                        start_n4,end_n4,alignLen_n4,teID_n4,uniqueTEID_n4,strand_n4 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i-4]
                        
                        # i -- current pos
                        start0,end0,alignLen0,teID0,uniqueTEID0,strand0 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i]
                        
                        # two downstream indexes
                        start1,end1,alignLen1,teID1,uniqueTEID1,strand1 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i+1]
                        start2,end2,alignLen2,teID2,uniqueTEID2,strand2 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i+2]
                        # comment out these two except for testing
                        start3,end3,alignLen3,teID3,uniqueTEID3,strand3 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i+3]
                        start4,end4,alignLen4,teID4,uniqueTEID4,strand4 = candidateSoloLTRCoordDict[ltrType][chrID][assemblyID][i+4]
                        # checking for teID matches that don't occur in two adjacent forward indexes
                        # print(teID_n2,start_n2,end_n2,teID_n1,start_n1,end_n1,teID0,start0,end0,teID1,start1,end1,teID2,start2,end2)
                        # if teID0 != teID1 and teID0 != teID2:
                        if teID0 != teID1 and teID0 != teID2 and teID0 != teID3 and teID0 != teID4:
                            # this comes before checking the upstream list slices because it considers the first position instead of starting with the 4th position in the list
                            if ltrType not in downstreamMatchDict:
                                downstreamMatchDict[ltrType] = {}
                            if chrID not in downstreamMatchDict[ltrType]:
                                downstreamMatchDict[ltrType][chrID] = {}
                            if assemblyID not in downstreamMatchDict[ltrType][chrID]:
                                downstreamMatchDict[ltrType][chrID][assemblyID] = {}
                            if uniqueTEID0 not in downstreamMatchDict[ltrType][chrID][assemblyID]:
                                downstreamMatchDict[ltrType][chrID][assemblyID][uniqueTEID0] = (start0,end0,alignLen0,teID0,strand0)
                            # not going to keep examples where there are downstream adjacent matches, so have this if-statement fall within earlier "downstream" if-statement
                            # this is checking if start and stop positions are less than current position. otherwise they correspond to the last entries in the last
                            # for the first two indexes of the list, searching i-1 and i-2 will look at the last two items of the list
                            # 
                            # if start_n1 < start0 and end_n1 < start0 and start_n2 < start0 and end_n2 < start0:
                            if start_n1 < start0 and end_n1 < start0 and start_n2 < start0 and end_n2 < start0 and start_n3 < start0 and end_n3 < start0 and start_n4 < start0 and end_n4 < start0:
                                # if teID0 != teID_n1 and teID0 != teID_n2:
                                if teID0 != teID_n1 and teID0 != teID_n2 and teID0 != teID_n3 and teID0 != teID_n4:
                                    # upstreamMatchDict
                                    if ltrType not in upstreamMatchDict:
                                        upstreamMatchDict[ltrType] = {}
                                    if chrID not in upstreamMatchDict[ltrType]:
                                        upstreamMatchDict[ltrType][chrID] = {}
                                    if assemblyID not in upstreamMatchDict[ltrType][chrID]:
                                        upstreamMatchDict[ltrType][chrID][assemblyID] = {}
                                    if uniqueTEID0 not in upstreamMatchDict[ltrType][chrID][assemblyID]:
                                        upstreamMatchDict[ltrType][chrID][assemblyID][uniqueTEID0] = (start0,end0,alignLen0,teID0,strand0)
    return(downstreamMatchDict,upstreamMatchDict)
                    

def createCountArray(scaffoldLen):
    countArray = [0]*scaffoldLen
    return(countArray)


def updateCountArray(countArray,start,stop):
    maxValue = 0
    # countArray[stop+1] += 1
    for i in range(start,stop):
        countArray[i] += 1
        if maxValue < countArray[i]:
            maxValue = countArray[i]
    return(countArray,maxValue)


def trackPositionInArray(countArray,start,stop,character):
    for i in range(start,stop):
        countArray[i] = character
    return(countArray)


def assessCandidateSoloLTRBoundariesWithIntactLTRs(chrLengthsDict,intactLTRCoordDict,candidateSoloLTRCoordDict,downstreamMatchDict,upstreamMatchDict,ltrType,soloWindowBoundary):
    filteredCandidateSoloLTRCoordDict = {}
    for assemblyID in chrLengthsDict:
        for chrID in chrLengthsDict[assemblyID]:
            chrLength   = chrLengthsDict[assemblyID][chrID]
            countArray  = createCountArray(chrLength)
            if ltrType in intactLTRCoordDict and chrID in intactLTRCoordDict[ltrType] and assemblyID in intactLTRCoordDict[ltrType][chrID]:
                for intactStart,intactStop,intactAlignLen,intactTEID,uniqueIntactTEID,intactStrand in intactLTRCoordDict[ltrType][chrID][assemblyID]:
                    if intactStart < intactStop:
                        countArray,maxValue = updateCountArray(countArray,intactStart,intactStop)
                    else:
                        print("unexpected coord orientation: ",ltrType,chrID,assemblyID,intactStart,intactStop)
                        sys.exit()
                        # intactLTRCoordDict['Copia'][chrID][assemblyID].append((int(start),int(end),alignLen,teID,uniqueTEID,strand))
                        # candidateSoloLTRCoordDict[ltrType][chrID][assemblyID].sort(key = lambda x:x[2], reverse=True)
            if ltrType in candidateSoloLTRCoordDict and chrID in candidateSoloLTRCoordDict[ltrType] and assemblyID in candidateSoloLTRCoordDict[ltrType][chrID]:
                for soloLTRStart,soloLTRStop,soloAlignLen,soloTEID,uniqueSoloTEID,soloStrand in candidateSoloLTRCoordDict[ltrType][chrID][assemblyID]:
                    if soloLTRStart < soloLTRStop:
                        countArray,maxValue = updateCountArray(countArray,soloLTRStart,soloLTRStop)
                    else:
                        print("unexpected coord orientation: ",ltrType,chrID,assemblyID,soloLTRStart,soloLTRStop)
                        sys.exit()
                    soloWindowStart = soloLTRStart - soloWindowBoundary
                    soloWindowStop = soloLTRStop + soloWindowBoundary
                    if soloWindowStart < 1:
                        soloWindowStart = 1
                    if soloWindowStop > chrLength:
                        soloWindowStop = chrLength
                    soloWindow5Prime = countArray[soloWindowStart:soloLTRStart]
                    soloWindow3Prime = countArray[soloLTRStop:soloWindowStop]
                    soloWindow5PrimeSet = set(soloWindow5Prime)
                    soloWindow3PrimeSet = set(soloWindow3Prime)
                    if len(soloWindow5PrimeSet) == 1:
                        if list(soloWindow5PrimeSet)[0] == 0:
                            if len(soloWindow3PrimeSet) == 1:
                                if list(soloWindow3PrimeSet)[0] == 0:
                                    if ltrType in downstreamMatchDict and ltrType in upstreamMatchDict:
                                        if chrID in downstreamMatchDict[ltrType] and chrID in upstreamMatchDict[ltrType]:
                                            if assemblyID in downstreamMatchDict[ltrType][chrID] and assemblyID in upstreamMatchDict[ltrType][chrID]:
                                                if uniqueSoloTEID in downstreamMatchDict[ltrType][chrID][assemblyID] and uniqueSoloTEID in upstreamMatchDict[ltrType][chrID][assemblyID]:
                                                    if ltrType not in filteredCandidateSoloLTRCoordDict:
                                                        filteredCandidateSoloLTRCoordDict[ltrType] = {}
                                                    if chrID not in filteredCandidateSoloLTRCoordDict[ltrType]:
                                                        filteredCandidateSoloLTRCoordDict[ltrType][chrID] = {}
                                                    if assemblyID not in filteredCandidateSoloLTRCoordDict[ltrType][chrID]:
                                                        filteredCandidateSoloLTRCoordDict[ltrType][chrID][assemblyID] = []
                                                    filteredCandidateSoloLTRCoordDict[ltrType][chrID][assemblyID].append((soloLTRStart,soloLTRStop,soloAlignLen,soloTEID,uniqueSoloTEID,soloStrand))
    return(filteredCandidateSoloLTRCoordDict)
                                                        

def assessCandidateSoloLTRBoundariesWithOtherSoloLTRs(chrLengthsDict,candidateSoloLTRCoordDict,downstreamMatchDict,upstreamMatchDict,ltrType,soloWindowBoundary):
    filteredCandidateSoloLTRCoordDict = {}
    # print("chrID,soloWindowBoundary,soloWindowStart,soloLTRStart,soloLTRStop,soloWindowStop")
    for assemblyID in chrLengthsDict:
        for chrID in chrLengthsDict[assemblyID]:
            chrLength   = chrLengthsDict[assemblyID][chrID]
            countArray  = createCountArray(chrLength)
            if ltrType in candidateSoloLTRCoordDict and chrID in candidateSoloLTRCoordDict[ltrType] and assemblyID in candidateSoloLTRCoordDict[ltrType][chrID]:
                # sort by length in descending order
                candidateSoloLTRCoordDict[ltrType][chrID][assemblyID].sort(key = lambda x:x[2], reverse=True)
                for soloLTRStart,soloLTRStop,soloAlignLen,soloTEID,uniqueSoloTEID,soloStrand in candidateSoloLTRCoordDict[ltrType][chrID][assemblyID]:
                    if soloLTRStart < soloLTRStop:
                        countArray,maxValue = updateCountArray(countArray,soloLTRStart,soloLTRStop)
                    else:
                        print("unexpected coord orientation: ",ltrType,chrID,assemblyID,intStart,intStop)
                        sys.exit()
                    soloWindowStart = soloLTRStart - soloWindowBoundary
                    soloWindowStop = soloLTRStop + soloWindowBoundary
                    if soloWindowStart < 1:
                        soloWindowStart = 1
                    if soloWindowStop > chrLength:
                        soloWindowStop = chrLength
                    soloWindow5Prime = countArray[soloWindowStart:soloLTRStart]
                    soloWindow3Prime = countArray[soloLTRStop:soloWindowStop]
                    soloWindow5PrimeSet = set(soloWindow5Prime)
                    soloWindow3PrimeSet = set(soloWindow3Prime)
                    #print(soloWindow5PrimeSet)
                    #print(soloWindow3PrimeSet)
                    # print(chrID,soloWindowBoundary,soloWindowStart,soloLTRStart,soloLTRStop,soloWindowStop)
                    ###
                    if len(soloWindow5PrimeSet) == 1:
                        if list(soloWindow5PrimeSet)[0] == 0:
                            if len(soloWindow3PrimeSet) == 1:
                                if list(soloWindow3PrimeSet)[0] == 0:
                                    # print(chrID,soloWindowBoundary,soloWindowStart,soloLTRStart,soloLTRStop,soloWindowStop)
                        
                                    ###
                                    if ltrType in downstreamMatchDict and ltrType in upstreamMatchDict:
                                        if chrID in downstreamMatchDict[ltrType] and chrID in upstreamMatchDict[ltrType]:
                                            if assemblyID in downstreamMatchDict[ltrType][chrID] and assemblyID in upstreamMatchDict[ltrType][chrID]:
                                                if uniqueSoloTEID in downstreamMatchDict[ltrType][chrID][assemblyID] and uniqueSoloTEID in upstreamMatchDict[ltrType][chrID][assemblyID]:
                                                    # candidateSoloLTRDataList.append((assemblyID,chrID,soloTEID,uniqueSoloTEID,soloLTRStart,soloLTRStop,soloAlignLen,soloStrand))
                                                    # filteredCandidateSoloLTRCoordDict
                                                    # for soloLTRStart,soloLTRStop,soloAlignLen,soloTEID,uniqueSoloTEID,soloStrand in candidateSoloLTRCoordDict[ltrType][chrID][assemblyID]:
                                                    if ltrType not in filteredCandidateSoloLTRCoordDict:
                                                        filteredCandidateSoloLTRCoordDict[ltrType] = {}
                                                    if chrID not in filteredCandidateSoloLTRCoordDict[ltrType]:
                                                        filteredCandidateSoloLTRCoordDict[ltrType][chrID] = {}
                                                    if assemblyID not in filteredCandidateSoloLTRCoordDict[ltrType][chrID]:
                                                        filteredCandidateSoloLTRCoordDict[ltrType][chrID][assemblyID] = []
                                                    filteredCandidateSoloLTRCoordDict[ltrType][chrID][assemblyID].append((soloLTRStart,soloLTRStop,soloAlignLen,soloTEID,uniqueSoloTEID,soloStrand))
                                                    
                    #else:
                    #    print(soloWindow5PrimeSet)
                    #    print(soloWindow3PrimeSet)
                    #    print(chrID,soloWindowBoundary,soloWindowStart,soloLTRStart,soloLTRStop,soloWindowStop) 
    return(filteredCandidateSoloLTRCoordDict)

                                                    
def assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,candidateSoloLTRCoordDict,internalFeatureCoordDict,downstreamMatchDict,upstreamMatchDict,ltrType,soloWindowBoundary):
    candidateSoloLTRLengthList = []
    candidateSoloLTRDataList = []
    featureCountDict = {}
    intChar = 1
    ltrChar = 2
    for assemblyID in chrLengthsDict:
        for chrID in chrLengthsDict[assemblyID]:
            chrLength   = chrLengthsDict[assemblyID][chrID]
            countArray  = createCountArray(chrLength)
            intCountArrayForRatio = createCountArray(chrLength)
            soloCountArrayForRatio = createCountArray(chrLength)
            if ltrType in internalFeatureCoordDict and chrID in internalFeatureCoordDict[ltrType] and assemblyID in internalFeatureCoordDict[ltrType][chrID]:
                for intStart,intStop,intAlignLen,intTEID,uniqueTEID,intStrand in internalFeatureCoordDict[ltrType][chrID][assemblyID]:
                    # checking coord orientation -- forward direction
                    if intStart < intStop:
                        countArray,maxValue = updateCountArray(countArray,intStart,intStop)
                        intCountArrayForRatio = trackPositionInArray(intCountArrayForRatio,intStart,intStop,intChar) 
                    else:
                        print("unexpected coord orientation: ",ltrType,chrID,assemblyID,intStart,intStop)
                        sys.exit()
            if ltrType in candidateSoloLTRCoordDict and chrID in candidateSoloLTRCoordDict[ltrType] and assemblyID in candidateSoloLTRCoordDict[ltrType][chrID]:
                for soloLTRStart,soloLTRStop,soloAlignLen,soloTEID,uniqueSoloTEID,soloStrand in candidateSoloLTRCoordDict[ltrType][chrID][assemblyID]:
                    # checking coord orientation
                    if soloLTRStart < soloLTRStop:
                        soloCountArrayForRatio = trackPositionInArray(soloCountArrayForRatio,soloLTRStart,soloLTRStop,ltrChar)
                        # this calculation finds the window flanking the solo-LTR, to check if there is internal sequence present
                        soloWindowStart = soloLTRStart - soloWindowBoundary
                        soloWindowStop = soloLTRStop + soloWindowBoundary
                        # if solo-LTR is located at edge of chromosome, it's possible that the window goes beyond length of chromosome
                        if soloWindowStart < 1:
                            soloWindowStart = 1
                        if soloWindowStop > chrLength:
                            soloWindowStop = chrLength
                        # python coords are not end-inclusive
                        # retrieve the section of the countArray the corresponds to the upstream or downstream LTR window
                        soloWindow5Prime = countArray[soloWindowStart:soloLTRStart]
                        soloWindow3Prime = countArray[soloLTRStop:soloWindowStop]
                        # get unique set
                        # x = [0,0,0,0,0,0,1,1]
                        # >>> set(x)
                        # {0, 1}
                        # this is reducing the countArraySubset into a set, which will be non-redundant characters, like a dict
                        soloWindow5PrimeSet = set(soloWindow5Prime)
                        soloWindow3PrimeSet = set(soloWindow3Prime)
                        # these if-statements are requiring that the set only contains 1 element and that element is zero (no other LTR elements present within 300 bp)
                        if len(soloWindow5PrimeSet) == 1:
                            if list(soloWindow5PrimeSet)[0] == 0:
                                if len(soloWindow3PrimeSet) == 1:
                                    if list(soloWindow3PrimeSet)[0] == 0:
                                        
                                        ### print(assemblyID,chrID,soloTEID,uniqueSoloTEID)
                                        ###
                                        ### require no matching downstream/upstream teIDs
                                        
                                        if ltrType in downstreamMatchDict and ltrType in upstreamMatchDict:
                                            if chrID in downstreamMatchDict[ltrType] and chrID in upstreamMatchDict[ltrType]:
                                                if assemblyID in downstreamMatchDict[ltrType][chrID] and assemblyID in upstreamMatchDict[ltrType][chrID]:
                                                    if uniqueSoloTEID in downstreamMatchDict[ltrType][chrID][assemblyID] and uniqueSoloTEID in upstreamMatchDict[ltrType][chrID][assemblyID]:
                                                        candidateSoloLTRLengthList.append(soloAlignLen)
                                                        candidateSoloLTRDataList.append((assemblyID,chrID,soloTEID,uniqueSoloTEID,soloLTRStart,soloLTRStop,soloAlignLen,soloStrand))
                    else:
                        print("unexpected coord orientation: ",ltrType,chrID,assemblyID,soloLTRStart,soloLTRStop,intStart,intStop)
                        sys.exit()
            # intChar = 1
            # ltrChar = 2
            intCount = intCountArrayForRatio.count(intChar)
            ltrCount = soloCountArrayForRatio.count(ltrChar)
            ltrIntRatio = float(ltrCount) / intCount
            if chrID not in featureCountDict:
                featureCountDict[chrID] = []
            featureCountDict[chrID].append(ltrIntRatio)
            # print(assemblyID,chrID,intCount,ltrCount,chrLength)
    return(candidateSoloLTRLengthList,candidateSoloLTRDataList)


def calculatePercentOverlap(start1,end1,start2,end2):
    shorterLength = min(end1-start1+1,end2-start2+1)
    overlapLength = 0
    # case I
    if start1 <= start2 and start2 <= end1:
        if end1 < end2:
            # case A    
            overlapLength = end1 - start2 + 1
        else:
            # case B
            overlapLength = end2 - start2 + 1
    # case II
    if start1 <= end2 and end2 <= end1:
        if start2 < start1:
            # case A
            overlapLength = end2 - start1 + 1
        else:
            # case B
            overlapLength = end2 - start2 + 1
    # case III
    if start2 <= start1 and start1 <= end2:
        if end2 < end1:
            # case A    
            overlapLength = end2 - start1 + 1
        else:
            # case B
            overlapLength = end1 - start1 + 1
    # case IV
    if start2 <= end1 and end1 <= end2:
        if start1 < start2:
            # case A
            overlapLength = end1 - start2 + 1
        else:
            # case B
            overlapLength = end1 - start1 + 1
    overlapFraction = float(overlapLength)/float(shorterLength)
    return(100.0*overlapFraction)


# candidateSoloLTRsWithDistances[ltrType][assemblyID][chrID][uniqueTEID1] = ('LTR',distanceBetweenFeatures)
def createCombinedStrandSpecificDict(plusDict,minusDict):
    combinedStrandSpecificDict = {}
    for ltrType in plusDict:
        ###
        if ltrType not in combinedStrandSpecificDict:
            combinedStrandSpecificDict[ltrType] = {}
        ###
        for assemblyID in plusDict[ltrType]:
            ###
            if assemblyID not in combinedStrandSpecificDict[ltrType]:
                combinedStrandSpecificDict[ltrType][assemblyID] = {}
            ###
            for chrID in plusDict[ltrType][assemblyID]:
                ###
                if chrID not in combinedStrandSpecificDict[ltrType][assemblyID]:
                    combinedStrandSpecificDict[ltrType][assemblyID][chrID] = {}
                ###
                for uniqueTEID in plusDict[ltrType][assemblyID][chrID]:
                    distanceFeatureType,distance = plusDict[ltrType][assemblyID][chrID][uniqueTEID]
                    if uniqueTEID not in combinedStrandSpecificDict[ltrType][assemblyID][chrID]:
                        combinedStrandSpecificDict[ltrType][assemblyID][chrID][uniqueTEID] = (distanceFeatureType,distance)
    for	ltrType	in minusDict:
        ###
        if ltrType not in combinedStrandSpecificDict:
            combinedStrandSpecificDict[ltrType] = {}
        ###
        for assemblyID in minusDict[ltrType]:
            ###
            if assemblyID not in combinedStrandSpecificDict[ltrType]:
                combinedStrandSpecificDict[ltrType][assemblyID] = {}
            ###
            for chrID in minusDict[ltrType][assemblyID]:
                ###
                if chrID not in combinedStrandSpecificDict[ltrType][assemblyID]:
                    combinedStrandSpecificDict[ltrType][assemblyID][chrID] = {}
                ###
                for uniqueTEID in minusDict[ltrType][assemblyID][chrID]:
                    distanceFeatureType,distance = minusDict[ltrType][assemblyID][chrID][uniqueTEID]
                    if uniqueTEID not in combinedStrandSpecificDict[ltrType][assemblyID][chrID]:
                        combinedStrandSpecificDict[ltrType][assemblyID][chrID][uniqueTEID] = (distanceFeatureType,distance)
    return(combinedStrandSpecificDict)
                

def writeCandidateSoloLTRs(candidateSoloLTRDataList,p,ltrType,strand,histBinWidth,specificAssemblyID,soloWindowBoundary):
    if ltrType  == 'Gypsy':
        newLTRType = 'Ty3'
    if ltrType == 'Copia':
        newLTRType = 'Ty1'
    OUT = open(specificAssemblyID + '_' + newLTRType + '_' + strand + '_binWidth' + str(histBinWidth) + '_soloLTRFlankingWindow' + str(soloWindowBoundary) + '_candidateSoloLTRData.tsv','w')
    OUT.write("# 95th percentile: %s\n" % (p))
    OUT.write("# Bin width: %s\n" % (histBinWidth))
    OUT.write("# Solo-LTR flanking window length: %s\n" % (soloWindowBoundary))
    # OUT.write("# Modal bin upper-bound: %s\n" % (maxBinSize))
    OUT.write("# assemblyID\tchrID\tltrType\tsoloTEID\tuniqueSoloTEID\tsoloLTRStart\tsoloLTRStop\talignLen\tstrand\n")
    for assemblyID,chrID,soloTEID,uniqueSoloTEID,soloLTRStart,soloLTRStop,alignLen,soloStrand in candidateSoloLTRDataList:
        if alignLen <= p:
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (assemblyID,chrID,newLTRType,soloTEID,uniqueSoloTEID,soloLTRStart,soloLTRStop,alignLen,soloStrand))

            

def createInternalSequenceLengthHist(internalSeqLengthDict,ltrType,strand,percentile,assemblyID):
    binWidth = 100
    if ltrType  == 'Gypsy':
        newLTRType = 'Ty3'
    if ltrType == 'Copia':
        newLTRType = 'Ty1'
    if ltrType in internalSeqLengthDict:
        for ltrType in internalSeqLengthDict:
            allLists = internalSeqLengthDict[ltrType]
            p = np.percentile(allLists, percentile)
            bins = np.arange(min(allLists),max(allLists),binWidth)
            fig, ax1 = plt.subplots(figsize=(6.25,5))
            plt.xlabel('Internal sequence lengths (kb)',size=16)
            plt.ylabel('Count',size=16)
            counts, bins, bars = plt.hist(allLists, bins = bins, histtype = 'step', color ='#0471A6')
            '''
            maxCount = 0
            for i in range(len(counts)):
                if counts[i] > maxCount:
                    maxCount = counts[i]
                    maxIndex = i
            avgMaxBin = float(bins[maxIndex] + bins[maxIndex+1]) / 2
            plt.axvline(avgMaxBin, linestyle='--', c='black')
            avgMaxBin = int(avgMaxBin)
            plt.text(avgMaxBin+binWidth, float(max(counts))/2, "average modal frequency, " + str(avgMaxBin) + " bp")
            '''
            maxCount = 0
            for i in range(len(counts)):
                if counts[i] > maxCount:
                    maxCount = counts[i]
                    maxIndex = i
            avgMaxBin = float(bins[maxIndex] + bins[maxIndex+1]) / 2
            plt.axvline(avgMaxBin, linestyle='--', c='black')
            avgMaxBin = int(avgMaxBin)
            #print("_internalSequenceLengthHist_binWidth: ",bins[maxIndex],bins[maxIndex+1],avgMaxBin)
            plt.text(avgMaxBin+binWidth, float(max(counts))/2, "average modal frequency, " + str(avgMaxBin) + " bp")

            plt.axvline(x = p, color = 'red', label = str(percentile) + "th percentile, " + str(p) + " bp")
            p = int(p)
            plt.annotate(str(percentile) + "th percentile, " + str(p) + " bp", (p,maxCount))
            
            x_ticks = ax1.get_xticks()/1000.0
            ax1.set_xticklabels(x_ticks.astype(int))
            plt.title(newLTRType)
            plt.tight_layout()
            plt.savefig(assemblyID + '_' + newLTRType + '_' + strand + '_internalSequenceLengthHist_binWidth' + str(binWidth) + '.png', dpi=600)
            plt.close()
    
            
def createFigure(candidateSoloLTRLengthList,percentile,p,ltrType,strand,assemblyID,soloWindowBoundary):
    if ltrType  == 'Gypsy':
        newLTRType = 'Ty3-LTR'
    if ltrType == 'Copia':
        newLTRType = 'Ty1-LTR'
    binWidth = 100
    allLists = candidateSoloLTRLengthList
    bins = np.arange(min(allLists),max(allLists),binWidth)
    fig, ax1 = plt.subplots(figsize=(7.5,6))
    plt.xlabel('Candidate solo-LTR lengths (kb)',size=12)
    plt.ylabel('Count',size=12)
    # plt.legend(frameon=False,fontsize=14)
    counts, bins, bars = plt.hist(candidateSoloLTRLengthList, bins = bins, histtype = 'step', color ='#0471A6')

    maxCount = 0
    for i in range(len(counts)):
        if counts[i] > maxCount:
            maxCount = counts[i]
            maxIndex = i
    avgMaxBin = float(bins[maxIndex] + bins[maxIndex+1]) / 2
    plt.axvline(avgMaxBin, linestyle='--', c='black')
    avgMaxBin = int(avgMaxBin)
    #print("candidateSoloLTRLengths: ",bins[maxIndex],bins[maxIndex+1],avgMaxBin)
    plt.text(avgMaxBin+binWidth, float(max(counts))/2, "average modal frequency, " + str(avgMaxBin) + " bp")
    
    plt.axvline(x = p, color = 'red', label = str(percentile) + "th percentile, " + str(p) + " bp")
    p = int(p)
    plt.annotate(str(percentile) + "th percentile, " + str(p) + " bp", (p,maxCount))
    
    # get ticks in units of 1 kb
    ax1.ticklabel_format(style='plain')
    x_ticks = ax1.get_xticks()/1000.0
    # y_ticks = ax.get_yticks()/1000000.0
    #print ticks
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_ticks.astype(float))
    # plt.tight_layout()
    plt.title(newLTRType)
    plt.savefig(assemblyID + '_' + newLTRType + '_' + strand + '_soloLTRFlankingWindow' + str(soloWindowBoundary) + '_candidateSoloLTRLengths.png', dpi=600)
    plt.close()
    
                    

usage = "Usage: " + sys.argv[0] + " <full TE anno gff file> <gene gff file> <distance histogram bin width> <assemblyID> <species or organism ID, i.e. cannabis or maize> <solo-LTR flanking window value, i.e. either 300 bp (LTR_retriever) or 5000 bp (LTR_MINER)>\n"
if len(sys.argv) != 7:
    print(usage)
    sys.exit()

teAnnoGFF = sys.argv[1]
geneGFF = sys.argv[2]
histBinWidth = sys.argv[3]
assemblyID = sys.argv[4]
organismID = sys.argv[5]
soloLTRFlankingWindow = sys.argv[6]

percentile = 95
histBinWidth = int(histBinWidth)
soloLTRFlankingWindow = int(soloLTRFlankingWindow)
soloWindowBoundary = soloLTRFlankingWindow

# chrLengthsDict = readGeneData(geneGFF)

# soloWindowBoundary = 300
# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2004-5-10-r79#MOESM5
# soloWindowBoundary = 5000

# chrLengthsDict = readGeneData(geneGFF)
if organismID == 'cannabis':
    chrLengthsDict = readCannabisGeneData(geneGFF)
elif organismID == 'maize':
    chrLengthsDict = readMaizeGeneData(geneGFF)
else:
    print("incorrect species/organism ID")
    sys.exit()

plus_candidateSoloLTRCoordDict,plus_internalFeatureCoordDict,plus_internalSeqLengthDict,plus_countDict,plus_intactLTRCoordDict = readTEAnnoGFF(teAnnoGFF,'+',assemblyID,organismID)
minus_candidateSoloLTRCoordDict,minus_internalFeatureCoordDict,minus_internalSeqLengthDict,minus_countDict,minus_intactLTRCoordDict = readTEAnnoGFF(teAnnoGFF,'-',assemblyID,organismID)

plus_downstreamMatchDict,plus_upstreamMatchDict = removeMatchingAdjacentEntries(plus_candidateSoloLTRCoordDict)
minus_downstreamMatchDict,minus_upstreamMatchDict = removeMatchingAdjacentEntries(minus_candidateSoloLTRCoordDict)

##### 
# Gypsy type, plus strand
plus_gypsyFilteredCandidateSoloLTRCoordDict1 = assessCandidateSoloLTRBoundariesWithIntactLTRs(chrLengthsDict,plus_intactLTRCoordDict,plus_candidateSoloLTRCoordDict,plus_downstreamMatchDict,plus_upstreamMatchDict,'Gypsy',soloWindowBoundary)
plus_gypsyFilteredCandidateSoloLTRCoordDict2 = assessCandidateSoloLTRBoundariesWithOtherSoloLTRs(chrLengthsDict,plus_gypsyFilteredCandidateSoloLTRCoordDict1,plus_downstreamMatchDict,plus_upstreamMatchDict,'Gypsy',soloWindowBoundary)
plus_gypsyCandidateSoloLTRLengthList,plus_gypsyCandidateSoloLTRDataList = assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,plus_gypsyFilteredCandidateSoloLTRCoordDict2,plus_internalFeatureCoordDict,plus_downstreamMatchDict,plus_upstreamMatchDict,'Gypsy',soloWindowBoundary)


# Gypsy type, minus strand
#minus_gypsyCandidateSoloLTRLengthList,minus_gypsyCandidateSoloLTRDataList = assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,minus_candidateSoloLTRCoordDict,minus_internalFeatureCoordDict,minus_downstreamMatchDict,minus_upstreamMatchDict,'Gypsy',soloWindowBoundary)
minus_gypsyFilteredCandidateSoloLTRCoordDict1 = assessCandidateSoloLTRBoundariesWithIntactLTRs(chrLengthsDict,minus_intactLTRCoordDict,minus_candidateSoloLTRCoordDict,minus_downstreamMatchDict,minus_upstreamMatchDict,'Gypsy',soloWindowBoundary)
minus_gypsyFilteredCandidateSoloLTRCoordDict2 = assessCandidateSoloLTRBoundariesWithOtherSoloLTRs(chrLengthsDict,minus_gypsyFilteredCandidateSoloLTRCoordDict1,minus_downstreamMatchDict,minus_upstreamMatchDict,'Gypsy',soloWindowBoundary)
minus_gypsyCandidateSoloLTRLengthList,minus_gypsyCandidateSoloLTRDataList = assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,minus_gypsyFilteredCandidateSoloLTRCoordDict2,minus_internalFeatureCoordDict,minus_downstreamMatchDict,minus_upstreamMatchDict,'Gypsy',soloWindowBoundary)

# combined Gypsy data
combinedGypsyLengthList = plus_gypsyCandidateSoloLTRLengthList + minus_gypsyCandidateSoloLTRLengthList
combinedGypsyDataList = plus_gypsyCandidateSoloLTRDataList + minus_gypsyCandidateSoloLTRDataList
gp = np.percentile(combinedGypsyLengthList, percentile)
#print("95th percentile solo-LTR (Gypsy, both strands) length (save values below this threshold): ",round(gp, 3))

#####
# Copia type, plus strand
# plus_copiaCandidateSoloLTRLengthList,plus_copiaCandidateSoloLTRDataList = assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,plus_candidateSoloLTRCoordDict,plus_internalFeatureCoordDict,plus_downstreamMatchDict,plus_upstreamMatchDict,'Copia',soloWindowBoundary)
plus_copiaFilteredCandidateSoloLTRCoordDict1 = assessCandidateSoloLTRBoundariesWithIntactLTRs(chrLengthsDict,plus_intactLTRCoordDict,plus_candidateSoloLTRCoordDict,plus_downstreamMatchDict,plus_upstreamMatchDict,'Copia',soloWindowBoundary)
plus_copiaFilteredCandidateSoloLTRCoordDict2 = assessCandidateSoloLTRBoundariesWithOtherSoloLTRs(chrLengthsDict,plus_copiaFilteredCandidateSoloLTRCoordDict1,plus_downstreamMatchDict,plus_upstreamMatchDict,'Copia',soloWindowBoundary)
plus_copiaCandidateSoloLTRLengthList,plus_copiaCandidateSoloLTRDataList = assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,plus_copiaFilteredCandidateSoloLTRCoordDict2,plus_internalFeatureCoordDict,plus_downstreamMatchDict,plus_upstreamMatchDict,'Copia',soloWindowBoundary)

# Copia type, minus strand 
# minus_copiaCandidateSoloLTRLengthList,minus_copiaCandidateSoloLTRDataList = assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,minus_candidateSoloLTRCoordDict,minus_internalFeatureCoordDict,minus_downstreamMatchDict,minus_upstreamMatchDict,'Copia',soloWindowBoundary)
minus_copiaFilteredCandidateSoloLTRCoordDict1 = assessCandidateSoloLTRBoundariesWithIntactLTRs(chrLengthsDict,minus_intactLTRCoordDict,minus_candidateSoloLTRCoordDict,minus_downstreamMatchDict,minus_upstreamMatchDict,'Copia',soloWindowBoundary)
minus_copiaFilteredCandidateSoloLTRCoordDict2 = assessCandidateSoloLTRBoundariesWithOtherSoloLTRs(chrLengthsDict,minus_copiaFilteredCandidateSoloLTRCoordDict1,minus_downstreamMatchDict,minus_upstreamMatchDict,'Copia',soloWindowBoundary)
minus_copiaCandidateSoloLTRLengthList,minus_copiaCandidateSoloLTRDataList = assessCandidateSoloLTRBoundariesWithInternalSequences(chrLengthsDict,minus_copiaFilteredCandidateSoloLTRCoordDict2,minus_internalFeatureCoordDict,minus_downstreamMatchDict,minus_upstreamMatchDict,'Copia',soloWindowBoundary)

# combined Copia data
combinedCopiaLengthList = plus_copiaCandidateSoloLTRLengthList + minus_copiaCandidateSoloLTRLengthList
combinedCopiaDataList = plus_copiaCandidateSoloLTRDataList + minus_copiaCandidateSoloLTRDataList
cp = np.percentile(combinedCopiaLengthList, percentile)
#print("95th percentile solo-LTR (Copia, both strands) length (save values below this threshold): ",round(cp, 3))

'''
# copia, plus
createInternalSequenceLengthHist(plus_internalSeqLengthDict,'Copia','plus',percentile,assemblyID)

# copia, minus
createInternalSequenceLengthHist(minus_internalSeqLengthDict,'Copia','minus',percentile,assemblyID)

# gypsy, plus
createInternalSequenceLengthHist(plus_internalSeqLengthDict,'Gypsy','plus',percentile,assemblyID)

# gypsy, minus
createInternalSequenceLengthHist(minus_internalSeqLengthDict,'Gypsy','minus',percentile,assemblyID)
'''

writeCandidateSoloLTRs(combinedGypsyDataList,gp,'Gypsy','bothStrands',histBinWidth,assemblyID,soloWindowBoundary)
createFigure(combinedGypsyLengthList,percentile,gp,'Gypsy','bothStrands',assemblyID,soloWindowBoundary)

writeCandidateSoloLTRs(combinedCopiaDataList,cp,'Copia','bothStrands',histBinWidth,assemblyID,soloWindowBoundary)
createFigure(combinedCopiaLengthList,percentile,cp,'Copia','bothStrands',assemblyID,soloWindowBoundary)

# plus_copiaCandidateSoloLTRDataList
# minus_copiaCandidateSoloLTRDataList
# plus_gypsyCandidateSoloLTRDataList
# minus_gypsyCandidateSoloLTRDataList
