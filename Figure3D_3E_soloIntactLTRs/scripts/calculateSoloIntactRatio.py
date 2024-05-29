import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


########
# MAIN #
########


def readFileList(fileList):
    dataList = []
    with open(fileList,'r') as F:
        for line in F:
            fileName = line.strip()
            dataList.append(fileName)
    return(dataList)


def readUnfilteredSoloLTRFile(unfilteredSoloLTRFile):
    soloLTRDict = {}
    with open(unfilteredSoloLTRFile,'r') as F:
        for line in F:
            if not line.startswith('assembly'):
                assemblyID,chrID,ltrType,soloTEID,uniqueSoloTEID,soloLTRStart,soloLTRStop,alignLen,score,percIden,strand = line.strip().split('\t')
                if ltrType not in soloLTRDict:
                    soloLTRDict[ltrType] = {}
                if assemblyID not in soloLTRDict[ltrType]:
                    soloLTRDict[ltrType][assemblyID] = {}
                if chrID not in soloLTRDict[ltrType][assemblyID]:
                    soloLTRDict[ltrType][assemblyID][chrID] = []
                soloLTRDict[ltrType][assemblyID][chrID].append((soloTEID,uniqueSoloTEID,int(soloLTRStart),int(soloLTRStop),int(alignLen)))
    return(soloLTRDict)


def readFilteredSoloLTRFile(filteredSoloLTRFile):
    soloLTRDict = {}
    with open(filteredSoloLTRFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                # assemblyID    chrID   ltrType soloTEID        uniqueSoloTEID  soloLTRStart    soloLTRStop     alignLen
                assemblyID,chrID,ltrType,soloTEID,uniqueSoloTEID,soloLTRStart,soloLTRStop,alignLen,strand = line.strip().split('\t')
                # print(assemblyID,chrID,ltrType,soloTEID,uniqueSoloTEID,soloLTRStart,soloLTRStop,alignLen)
                if ltrType not in soloLTRDict:
                    soloLTRDict[ltrType] = {}
                if assemblyID not in soloLTRDict[ltrType]:
                    soloLTRDict[ltrType][assemblyID] = {}
                if chrID not in soloLTRDict[ltrType][assemblyID]:
                    soloLTRDict[ltrType][assemblyID][chrID] = []
                soloLTRDict[ltrType][assemblyID][chrID].append((soloTEID,uniqueSoloTEID,int(soloLTRStart),int(soloLTRStop),int(alignLen)))
    return(soloLTRDict)


def readTransposonGFF(transposonDataList,organismID):
    gffData = {}
    gffCoordDict = {}
    minAlignmentLength = 100
    minPercentIdentity = 0.8
    for gffFile in transposonDataList:
        with open(gffFile,'r') as F:
            for line in F:
                if not line.startswith('#'):
                    scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                    start = int(start)
                    end = int(end)
                    seqLen = end - start + 1
                    if organismID == 'cannabis':
                        assemblyID,chrID = scaffoldID.split('.')
                    elif organismID == 'maize':
                        basename = os.path.basename(gffFile)
                        assemblyID,suffixID = basename.strip().split('-REFERENCE')
                        chrID = scaffoldID
                    if 'chr' in scaffoldID:
                        if 'LTR_retrotransposon' in feature:
                            if 'Method=structural' in attribute:
                                # ['ID=LTRRT_44', 'Parent=repeat_region_44', 'Name=TE_00003284', 'Classification=LTR/Copia', 'Sequence_ontology=SO:0002264', 'ltr_identity=0.9802', 'Method=structural', 'motif=TGCA', 'tsd=TTAAA']
                                newAttribute = attribute.split(';')
                                suffix,ltrIdentity = newAttribute[5].split('ltr_identity=')
                                ltrIdentity = ltrIdentity.strip()
                                ltrIdentity = float(ltrIdentity)
                                ltrDist = 1-ltrIdentity
                                if ltrIdentity >= minPercentIdentity:
                                    if seqLen >= minAlignmentLength:
                                        # print(scaffoldID,source,feature,start,end,score,strand,frame,attribute)
                                        if 'Gypsy' in feature:
                                            if 'Ty3' not in gffData:
                                                gffData['Ty3'] = {}
                                            if chrID not in gffData['Ty3']:
                                                gffData['Ty3'][chrID] = {}
                                            if assemblyID not in gffData['Ty3'][chrID]:
                                                gffData['Ty3'][chrID][assemblyID] = []
                                            gffData['Ty3'][chrID][assemblyID].append(ltrDist)

                                            if 'Ty3' not in gffCoordDict:
                                                gffCoordDict['Ty3'] = {}
                                            if chrID not in gffCoordDict['Ty3']:
                                                gffCoordDict['Ty3'][chrID] = {}
                                            if assemblyID not in gffCoordDict['Ty3'][chrID]:
                                                gffCoordDict['Ty3'][chrID][assemblyID] = []
                                            gffCoordDict['Ty3'][chrID][assemblyID].append((int(start),int(end),seqLen))
                                        elif 'Copia' in feature:
                                            if 'Ty1' not in gffData:
                                                gffData['Ty1'] = {}
                                            if chrID not in gffData['Ty1']:
                                                gffData['Ty1'][chrID] = {}
                                            if assemblyID not in gffData['Ty1'][chrID]:
                                                gffData['Ty1'][chrID][assemblyID] = []
                                            gffData['Ty1'][chrID][assemblyID].append(ltrDist)
                                
                                            if 'Ty1' not in gffCoordDict:
                                                gffCoordDict['Ty1'] = {}
                                            if chrID not in gffCoordDict['Ty1']:
                                                gffCoordDict['Ty1'][chrID] = {}
                                            if assemblyID not in gffCoordDict['Ty1'][chrID]:
                                                gffCoordDict['Ty1'][chrID][assemblyID] = []
                                            gffCoordDict['Ty1'][chrID][assemblyID].append((int(start),int(end),seqLen))
    return(gffData,gffCoordDict)


def calculateSoloIntactRatio(gffData,soloLTRDict,specificLTRType):
    dataDict = {}
    for ltrType in gffData:
        if ltrType == specificLTRType:
            OUT = open(ltrType + '_soloIntactRatio.tsv','w')
            OUT.write("ltrType\tassemblyID\tchrID\tsoloCount\tintactCount\tSIRatio\n")
            for chrID in gffData[ltrType]:
                for assemblyID in gffData[ltrType][chrID]:
                    intactCount = len(gffData[ltrType][chrID][assemblyID])
                    if ltrType in soloLTRDict:
                        soloCount = len(soloLTRDict[ltrType][assemblyID][chrID])
                        SIRatio = float(soloCount) / intactCount
                        # print(ltrType,assemblyID,chrID,soloCount,intactCount,SIRatio)
                        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ltrType,assemblyID,chrID,soloCount,intactCount,SIRatio))
                        if chrID not in dataDict:
                            dataDict[chrID] = []
                        dataDict[chrID].append(SIRatio)
    return(dataDict)

                    
def getStdDev(dataList):
    dataStdDev = round(np.std(dataList), 2)
    return(dataStdDev)


def createDataArray(specificDataDict):
    fullDataArray = []
    fullStdDevList = []
    idList = []
    for infoID in specificDataDict:
        idList.append(infoID)
        dataList = specificDataDict[infoID]
        dataStdDev = getStdDev(dataList)
        fullDataArray.append(np.asarray(dataList))
        fullStdDevList.append(dataStdDev)
    # https://stackoverflow.com/questions/19931975/sort-multiple-lists-simultaneously
    idList,fullDataArray,fullStdDevList = map(list, zip(*sorted(zip(idList, fullDataArray, fullStdDevList), reverse=False)))
    return(fullDataArray,fullStdDevList,idList)


# https://github.com/oushujun/EDTA/issues/233
# "You may grep LTR_retrotransposon *.fasta.mod.EDTA.intact.gff3|less -S to get only intact LTRs."
# "For intact LTR information in EDTA results, you can find the identity info from the last column of the gff3 file. You can use this to calculate the age yourself, with T = K / 2µ = (1 - identity) / 2µ, where µ is the mutation rate of your species in the unit of per bp per year. By default, LTR_retriever uses µ = 1.3e-8 from rice."
def createBoxplot(dataArray,stdDevList,idList,ltrType,infoLabel,soloFlankingWindow):
    # Create a figure instance for subs
    plt.rcParams["figure.figsize"] = [4,2]
    #colors = ['#278B9A','#278B9A','#E75B64','#E75B64','#D8AF39','#D8AF39','#33658A','#33658A']
    colors = ['#278B9A']*len(stdDevList)
    
    # Create the boxplot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    #labelList = stdDevList
    labelList = []
    exonCount = 0
    for infoID in idList:
        label = infoID
        labelList.append(label)

    bp = ax.boxplot(dataArray, labels = labelList, patch_artist=True, showfliers=False, zorder=0, showmeans=True, meanline=True)
    ax.tick_params(axis='x', labelsize=6)
    plt.xticks(rotation=45, ha="right")

    plt.ylabel('Solo:intact LTR ratio \n (all assemblies)',size=10)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(0.6)

    for i in range(len(dataArray)):
        y1 = dataArray[i]
        x1 = np.random.normal(i+1, 0.02, len(y1))
        plt.plot(x1, y1, 'black', alpha=0.3, marker="o", markersize=3, linestyle=' ',
                 fillstyle='full',
                 markeredgecolor='black',
                 markeredgewidth=0.0)

    plt.title(ltrType)
    plt.savefig(ltrType + '_' + infoLabel + '_soloFlankingWindow' + str(soloFlankingWindow) + '_soloIntactRatioBoxplot.png', bbox_inches='tight', dpi=600)
    plt.savefig(ltrType + '_' + infoLabel + '_soloFlankingWindow' + str(soloFlankingWindow) + '_soloIntactRatioBoxplot.pdf', bbox_inches='tight')
    plt.savefig(ltrType + '_' + infoLabel + '_soloFlankingWindow' + str(soloFlankingWindow) + '_soloIntactRatioBoxplot.svg', bbox_inches='tight')
    plt.close()



usage = "Usage: " + sys.argv[0] + " <unfiltered solo ltr file> <filtered solo ltr file> <intact gff file list> <solo flanking window> <organism/species ID>\n"
if len(sys.argv) != 6:
    print(usage)
    sys.exit()

unfilteredSoloLTRFile = sys.argv[1]
filteredSoloLTRFile = sys.argv[2]
transposonGFFFileList = sys.argv[3]
soloFlankingWindow = sys.argv[4]
organismID = sys.argv[5]

unfilteredSoloLTRDict = readUnfilteredSoloLTRFile(unfilteredSoloLTRFile)
filteredSoloLTRDict = readFilteredSoloLTRFile(filteredSoloLTRFile)

transposonDataList = readFileList(transposonGFFFileList)
gffData,gffCoordDict = readTransposonGFF(transposonDataList,organismID)

unfilteredDataDictTy3 = calculateSoloIntactRatio(gffData,unfilteredSoloLTRDict,'Ty3')
unfilteredDataDictTy1 = calculateSoloIntactRatio(gffData,unfilteredSoloLTRDict,'Ty1')
filteredDataDictTy3 = calculateSoloIntactRatio(gffData,filteredSoloLTRDict,'Ty3')
filteredDataDictTy1 = calculateSoloIntactRatio(gffData,filteredSoloLTRDict,'Ty1')

# uf = unfiltered
# f = filtered
uf_dataArrayTy3,uf_stdDevListTy3,uf_idListTy3 = createDataArray(unfilteredDataDictTy3)
uf_dataArrayTy1,uf_stdDevListTy1,uf_idListTy1 = createDataArray(unfilteredDataDictTy1)

f_dataArrayTy3,f_stdDevListTy3,f_idListTy3 = createDataArray(filteredDataDictTy3)
f_dataArrayTy1,f_stdDevListTy1,f_idListTy1 = createDataArray(filteredDataDictTy1)

createBoxplot(uf_dataArrayTy3,uf_stdDevListTy3,uf_idListTy3,'Ty3','unfiltered','None')
createBoxplot(uf_dataArrayTy1,uf_stdDevListTy1,uf_idListTy1,'Ty1','unfiltered','None')

createBoxplot(f_dataArrayTy3,f_stdDevListTy3,f_idListTy3,'Ty3','filtered',soloFlankingWindow)
createBoxplot(f_dataArrayTy1,f_stdDevListTy1,f_idListTy1,'Ty1','filtered',soloFlankingWindow)
