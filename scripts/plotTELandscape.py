import sys, re, os, math, statistics
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from pandas import DataFrame
from mpl_toolkits.axes_grid1 import make_axes_locatable


########
# MAIN #
########


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


def readOtherAnnotationBedFiles(otherAnnotationBedFile,moleculeID,orientationDict,assemblyLengthDict):
    annotationCoordDict = {}
    with open(otherAnnotationBedFile,'r') as F:
        for line in F:
            if '#' not in line:
                fullChromID,start,stop,annotationID,fifthCol,strand = line.strip().split('\t')
                start = int(start)
                stop = int(stop)
                assemblyID,chrID = fullChromID.split('.')
                #SODLb.chr4      71136148        71137277        SODLb_0 0       -
                #SODLb.chr4      71137302        71140750        SODLb_1 0       -
                if assemblyID in assemblyLengthDict and chrID in assemblyLengthDict[assemblyID]:
                    #print(fullChromID,start,stop,annotationID)
                    chromLength = assemblyLengthDict[assemblyID][chrID]
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
                        #print(start,stop)
                    else:
                        print(assemblyID,chrID," not in orientationDict")
                        sys.exit()
                if fullChromID not in annotationCoordDict:
                    annotationCoordDict[fullChromID] = {}
                if moleculeID not in annotationCoordDict[fullChromID]:
                    annotationCoordDict[fullChromID][moleculeID] = []
                annotationCoordDict[fullChromID][moleculeID].append((int(start),int(stop),annotationID))
    return(annotationCoordDict)


def readBedFile(synthaseBedFile,orientationDict,assemblyLengthDict):
    synthaseCoordDict = {}
    with open(synthaseBedFile,'r') as F:
        for line in F:
            if '#' not in line:
                fullChromID,start,stop,synthaseID,fifthCol,strand = line.strip().split('\t')
                start = int(start)
                stop = int(stop)
                assemblyID,chrID = fullChromID.split('.')
                items = synthaseID.split('_')
                moleculeID = items[1]
                if assemblyID in assemblyLengthDict and chrID in assemblyLengthDict[assemblyID]:
                    chromLength = assemblyLengthDict[assemblyID][chrID]
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
                        print(assemblyID,chrID," not in orientationDict")
                        sys.exit()
                if fullChromID not in synthaseCoordDict:
                    synthaseCoordDict[fullChromID] = {}
                if moleculeID not in synthaseCoordDict[fullChromID]:
                    synthaseCoordDict[fullChromID][moleculeID] = []
                synthaseCoordDict[fullChromID][moleculeID].append((start,stop,synthaseID))
    return(synthaseCoordDict)


def readSoloIntactRatioFile(soloIntactRatioFile,windowStart,windowStop):
    specificChromDataDict = {}
    with open(soloIntactRatioFile,'r') as F:
        for line in F:
            if 'ltrType' not in line:
                ltrType,assemblyID,chrID,genomicPosition,soloCount,intactCount,siRatio = line.strip().split('\t')
                # Ty1     EH23a   chr1    0.0     1       4       0.25
                genomicPosition = float(genomicPosition)
                fullChromID = assemblyID + '.' + chrID
                genomicPosition = int(float(genomicPosition))
                soloCount = int(soloCount)
                intactCount = int(intactCount)
                if genomicPosition >= windowStart and genomicPosition < windowStop:
                    if siRatio != 'nan':
                        siRatio = float(siRatio)
                    if assemblyID not in specificChromDataDict:
                        specificChromDataDict[assemblyID] = {}
                    if chrID not in specificChromDataDict[assemblyID]:
                        specificChromDataDict[assemblyID][chrID] = {}
                    if ltrType not in specificChromDataDict[assemblyID][chrID]:
                        specificChromDataDict[assemblyID][chrID][ltrType] = []
                    specificChromDataDict[assemblyID][chrID][ltrType].append((genomicPosition,soloCount,intactCount,siRatio))
    for assemblyID in specificChromDataDict:
        for chrID in specificChromDataDict[assemblyID]:
            for ltrType in specificChromDataDict[assemblyID][chrID]:
                specificChromDataDict[assemblyID][chrID][ltrType].sort(key=lambda x:x[0], reverse=False)
    return(specificChromDataDict)


def readChromosomeLengthData(chromosomeLengths):
    chromosomeLengthDict = {}
    with open(chromosomeLengths,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,chromLength = line.strip().split('\t')
                if assemblyID not in chromosomeLengthDict:
                    chromosomeLengthDict[assemblyID] = {}
                if chromID not in chromosomeLengthDict[assemblyID]:
                    chromosomeLengthDict[assemblyID][chromID] = int(chromLength)
    return(chromosomeLengthDict)


def calculateLastWindowLength(chromosomeLengthDict,windowSize):
    lastWindowRemainderDict = {}
    for assemblyID in chromosomeLengthDict:
        for chromID in chromosomeLengthDict[assemblyID]:
            chromLength = chromosomeLengthDict[assemblyID][chromID]
            lastWindowLength = math.fmod(chromLength, windowSize)
            if assemblyID not in lastWindowRemainderDict:
                lastWindowRemainderDict[assemblyID] = {}
            if chromID not in lastWindowRemainderDict[assemblyID]:
                lastWindowRemainderDict[assemblyID][chromID] = lastWindowLength
    return(lastWindowRemainderDict)


def readGeneData(geneData,windowStart,windowStop):
    geneDensity = {}
    with open(geneData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,position,geneCount = line.strip().split('\t')
                position = float(position)
                if position >= windowStart and position < windowStop:
                    if assemblyID not in geneDensity:
                        geneDensity[assemblyID] = {}
                    if chromID not in geneDensity[assemblyID]:
                        geneDensity[assemblyID][chromID] = []
                    geneDensity[assemblyID][chromID].append((float(position),int(geneCount)))
    for assemblyID in geneDensity:
        for chromID in geneDensity[assemblyID]:
            geneDensity[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
    return(geneDensity)


def readTETimeData(timeData,windowStart,windowStop):
    timeDataDict = {}
    with open(timeData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,featureID,timeID,chromStartPosition,specificTECount,teCount,totalTELenPerWindow,percentWindow,averageTime,specificTimeTEPercent = line.strip().split('\t')
                chromStartPosition = float(chromStartPosition)
                if chromStartPosition >= windowStart and chromStartPosition < windowStop:
                    if assemblyID not in timeDataDict:
                        timeDataDict[assemblyID] = {}
                    if chromID not in timeDataDict[assemblyID]:
                        timeDataDict[assemblyID][chromID] = {}
                    if featureID not in timeDataDict[assemblyID][chromID]:
                        timeDataDict[assemblyID][chromID][featureID] = {}
                    if timeID not in timeDataDict[assemblyID][chromID][featureID]:
                        timeDataDict[assemblyID][chromID][featureID][timeID] = []
                    timeDataDict[assemblyID][chromID][featureID][timeID].append((float(chromStartPosition),int(teCount),int(totalTELenPerWindow),float(percentWindow),float(averageTime),float(specificTimeTEPercent)))
    for assemblyID in timeDataDict:
        for chromID in timeDataDict[assemblyID]:
            for featureID in timeDataDict[assemblyID][chromID]:
                for timeID in timeDataDict[assemblyID][chromID][featureID]:
                    timeDataDict[assemblyID][chromID][featureID][timeID].sort(key=lambda x:x[0], reverse=False)
    return(timeDataDict)
                

def readTEDensityData(teData,windowStart,windowStop):
    densityData = {}
    featureIDs = {}
    with open(teData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,featureID,position,teCount,totalTELenPerWindow,percentWindow = line.strip().split('\t')
                position = float(position)
                if position >= windowStart and position < windowStop:
                    if featureID not in featureIDs:
                        featureIDs[featureID] = 1
                    if assemblyID not in densityData:
                        densityData[assemblyID] = {}
                    if chromID not in densityData[assemblyID]:
                        densityData[assemblyID][chromID] = {}
                    if featureID not in densityData[assemblyID][chromID]:
                        densityData[assemblyID][chromID][featureID] = []
                    densityData[assemblyID][chromID][featureID].append((float(position),int(teCount),int(totalTELenPerWindow),float(percentWindow)))
    for assemblyID in densityData:
        for chromID in densityData[assemblyID]:
            for featureID in densityData[assemblyID][chromID]:
                densityData[assemblyID][chromID][featureID].sort(key=lambda x:x[0], reverse=False)
    return(densityData,featureIDs)


def readMethylationData(methylationData,windowStart,windowStop):
    methylationDensity = {}
    with open(methylationData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,position,averageProbability = line.strip().split('\t')
                position = float(position)
                if position >= windowStart and position < windowStop:
                    if assemblyID not in methylationDensity:
                        methylationDensity[assemblyID] = {}
                    if chromID not in methylationDensity[assemblyID]:
                        methylationDensity[assemblyID][chromID] = []
                    methylationDensity[assemblyID][chromID].append((float(position),float(averageProbability)))
                    #print('#',position,windowStart,windowStop,assemblyID,chromID)
    for assemblyID in methylationDensity:
        for chromID in methylationDensity[assemblyID]:
            methylationDensity[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
    return(methylationDensity)


def readGCData(gcData,windowStart,windowStop):
    gcDataDict = {}
    with open(gcData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,position,observed_vs_expected,gc_content = line.strip().split('\t')
                position = float(position)
                gc_content = float(gc_content)
                gc_content = gc_content*100
                if position >= windowStart and position < windowStop:
                    if assemblyID not in gcDataDict:
                        gcDataDict[assemblyID] = {}
                    if chromID not in gcDataDict[assemblyID]:
                        gcDataDict[assemblyID][chromID] = []
                    gcDataDict[assemblyID][chromID].append((float(position),float(observed_vs_expected),gc_content))
                    #print('#',position,windowStart,windowStop,assemblyID,chromID)
    for assemblyID in gcDataDict:
        for chromID in gcDataDict[assemblyID]:
            gcDataDict[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
    return(gcDataDict)


def readVariantData(variantData,windowStart,windowStop):
    variantDensityDict = {}
    with open(variantData,'r') as F:
        for line in F:
            if 'assemblyID' not in line:
                assemblyID,chromID,position,variantCount,totalVariantLength,totalVariantPercent = line.strip().split('\t')
                position = float(position)
                if position >= windowStart and position < windowStop:
                    if assemblyID not in variantDensityDict:
                        variantDensityDict[assemblyID] = {}
                    if chromID not in variantDensityDict[assemblyID]:
                        variantDensityDict[assemblyID][chromID] = []
                    variantDensityDict[assemblyID][chromID].append((float(position),int(variantCount),int(totalVariantLength),float(totalVariantPercent)))
    for assemblyID in variantDensityDict:
        for chromID in variantDensityDict[assemblyID]:
            variantDensityDict[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
    return(variantDensityDict)



def plotCondensed(densityData,compiledTimes,methylationDensity,chromosomeLengthDict,gcDataDict,geneDensity,lastWindowRemainderDict,featureType,windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict,ty1_specificChromDataDict,ty3_specificChromDataDict,teColorMap,ltrColorMap,windowStart,windowStop,synthaseCoordDict,specificTE):
    for assemblyID in densityData:
        for chromID in densityData[assemblyID]:
            fullChromID = assemblyID + '.' + chromID
            seqLen = chromosomeLengthDict[assemblyID][chromID]
            lastWindowLength = lastWindowRemainderDict[assemblyID][chromID]
            fig, ax1 = plt.subplots(figsize=(25,6.25))
            #fig, (ax1, ax6, ax7, ax8, ax9, ax10, ax13) = plt.subplots(nrows=7, sharex=True, figsize=(25,12.5), gridspec_kw=dict(height_ratios=[3, 1.5, 1.5, 1, 1, 1, 1]))
            #fig, (ax1,ax7) = plt.subplots(figsize=(25,12.5), nrows=2, sharex=True)
            plt.rcParams['font.size'] = 16
            legendDict1 = {}
            for featureID in densityData[assemblyID][chromID]:
                featureArray = np.array(densityData[assemblyID][chromID][featureID])
                if featureType == 'LTRs':
                    colorMap = ltrColorMap
                else:
                    colorMap = teColorMap
                for i in range(len(featureArray)-1):
                    position,teCount,teTotalLen,totalTEPercent = featureArray[i]
                    ax1.bar(position,totalTEPercent,windowSize, align='edge', label=featureID, alpha=1, edgecolor=colorMap[featureID], color='None', linewidth=4)
                lastPosition,lastTECount,lastTotalTELen,lastTotalTEPercent = featureArray[-1]
                testStop = lastPosition+windowSize
                if testStop > seqLen:
                    ax1.bar(lastPosition,lastTotalTEPercent,lastWindowLength, align='edge', label=featureID, alpha=1, edgecolor=colorMap[featureID], color='None', linewidth=4)
                else:
                    ax1.bar(lastPosition,lastTotalTEPercent,windowSize, align='edge', label=featureID, alpha=1, edgecolor=colorMap[featureID], color='None', linewidth=4)
                ax1.spines['left'].set_color('black')
                ax1.spines['right'].set_color('black')
                ax1.tick_params(axis='y', colors='black')
                ax1.set_ylim(0,100)
                ax1.set_ylabel('Transposon percent (%)')
                ax1.margins(0)
                legendInfo = mlines.Line2D([], [], color=colorMap[featureID], marker='s',linestyle="None",markersize=10, label=featureID, alpha=1.0, fillstyle='full', markeredgecolor=colorMap[featureID], markeredgewidth=0.0)
                ticks = ax1.get_xticks()/1000.0
                ax1.set_xticks(ticks)
                ax1.set_xticklabels(ticks.astype(int))
                if featureID not in legendDict1:
                    legendDict1[featureID] = legendInfo
            synthaseColors = {'THCAS':'blue', 'CBDAS':'red', 'CBCAS':'orange'}
            if fullChromID in synthaseCoordDict:
                for moleculeType in synthaseCoordDict[fullChromID]:
                    for synthaseStart,synthaseStop,synthaseID in synthaseCoordDict[fullChromID][moleculeType]:
                        if synthaseStart >= windowStart and synthaseStart < windowStop:
                            ax1.vlines(synthaseStart, 0, 100, linewidth=4, color=synthaseColors[moleculeType])
                        else:
                            ax1.axvspan(synthaseStart,synthaseStop, ymin=0, ymax=100, color=synthaseColors[moleculeType], alpha=1, linestyle='solid', linewidth=2)
                        if len(synthaseCoordDict[fullChromID][moleculeType]) > 1:
                            synthaseLegendInfo = mlines.Line2D([], [], color=synthaseColors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + ', ' + str(len(synthaseCoordDict[fullChromID][moleculeType])) + ' copies', alpha=1.0, fillstyle='full', markeredgecolor=synthaseColors[moleculeType], markeredgewidth=0.0, linewidth=4)
                        else:
                            synthaseLegendInfo = mlines.Line2D([], [], color=synthaseColors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + '\n' + str(synthaseStart) + '-' + str(synthaseStop) + ' bp', alpha=1.0, fillstyle='full', markeredgecolor=synthaseColors[moleculeType], markeredgewidth=0.0, linewidth=4)
                        if moleculeType not in legendDict1:
                            legendDict1[moleculeType] = synthaseLegendInfo
                            
            ax2 = ax1.twinx()
            # if position >= windowStart and position < windowStop:
            if assemblyID in methylationDensity and chromID in methylationDensity[assemblyID]:
                methylationArray = np.array(methylationDensity[assemblyID][chromID])
                for j in range(len(methylationArray)-2):
                    position,methylationProb = methylationArray[j]
                    position = position + (windowSize/2)
                    nextPosition,nextMethylationProb = methylationArray[j+1]
                    nextPosition = nextPosition + (windowSize/2)
                    ax2.plot((position,nextPosition),(methylationProb,nextMethylationProb), color='black', marker='o', linewidth=4, markersize=6)
                secondToLastPosition,secondToLastMethylationProb = methylationArray[-2]
                lastPos,lastMethylationProb = methylationArray[-1]
                secondToLastPosition = secondToLastPosition + (windowSize/2)
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPosition = lastPos + (lastWindowLength/2)
                else:
                    lastPosition = lastPos + (windowSize/2)
                ax2.plot((secondToLastPosition,lastPosition),(secondToLastMethylationProb,lastMethylationProb), color='black', marker='o', linewidth=4, markersize=6)
                ax2.set_ylabel('Average methylation probability')
                ax2.margins(0)
                ax2.set_ylim(0,100)
                ticks = ax2.get_xticks()/1000.0
                ax2.set_xticks(ticks)
                ax2.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='black', marker='o',markersize=6, label='Methylation probability', alpha=1.0, fillstyle='full', markeredgecolor='black', markeredgewidth=0.0)
                if 'Methylation_probability' not in legendDict1:
                    legendDict1['Methylation_probability'] = legendInfo
            ax3 = ax1.twinx()
            ax3.grid(False)
            ax3.spines.right.set_position(("axes", 1.08))
            if assemblyID in gcDataDict and chromID in gcDataDict[assemblyID]:
                positions,observed_expected,gcContent = zip(*gcDataDict[assemblyID][chromID])
                smallestPos = min(positions)
                largestPos = max(positions)
                largestObservedExpected = max(observed_expected)
                positionArray = np.array(positions)
                obs_expArray = np.array(observed_expected)
                for k in range(len(obs_expArray)-2):
                    position = positionArray[k]
                    position = position + (windowSize/2)
                    obs_exp = obs_expArray[k]
                    nextPosition = positionArray[k+1]
                    nextPosition = nextPosition + (windowSize/2)
                    next_obs_exp = obs_expArray[k+1]
                    p3, = ax3.plot((position,nextPosition),(obs_exp,next_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=6)
                secondToLastPos = positionArray[-2]
                lastPos = positionArray[-1]
                secondToLastPos = secondToLastPos + (windowSize/2)
                secondToLast_obs_exp = obs_expArray[-2]
                last_obs_exp = obs_expArray[-1]
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPos = lastPos + (lastWindowLength/2)
                    ax3.hlines(y=obs_exp_threshold, xmin=smallestPos, xmax=largestPos+lastWindowLength, linewidth=3, color='#E0A458')
                else:
                    lastPos = lastPos + (windowSize/2)
                    ax3.hlines(y=obs_exp_threshold, xmin=smallestPos, xmax=largestPos+windowSize, linewidth=3, color='#E0A458')
                ax3.plot((secondToLastPos,lastPos),(secondToLast_obs_exp,last_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=6)
                ax3.spines['right'].set_color('#E0A458')
                ax3.tick_params(axis='y', colors='#E0A458')
                ax3.set_ylabel('Observed/expected CpG')
                ax3.yaxis.get_label().set_color(p3.get_color())
                ax3.margins(0)
                if largestObservedExpected > 1:
                    ax3.set_ylim(0,)
                else:
                    ax3.set_ylim(0,1)
                ticks = ax3.get_xticks()/1000.0
                ax3.set_xticks(ticks)
                ax3.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='#E0A458', marker='o',markersize=6, label='Observed/expected CpG', alpha=1.0, fillstyle='full', markeredgecolor='#E0A458', markeredgewidth=0.0)
                if 'Observed/expected CpG' not in legendDict1:
                    legendDict1['Observed/expected CpG'] = legendInfo
            ax4 = ax1.twinx()
            ax4.grid(False)
            ax4.spines.left.set_position(("axes", -0.08))
            if assemblyID in gcDataDict and chromID in gcDataDict[assemblyID]:
                positions,observed_expected,gcContent = zip(*gcDataDict[assemblyID][chromID])
                positionArray = np.array(positions)
                gcContentArray = np.array(gcContent)
                for l in range(len(gcContentArray)-2):
                    position = positionArray[l]
                    position = position + (windowSize/2)
                    nextPosition = positionArray[l+1]
                    nextPosition = nextPosition + (windowSize/2)
                    gc = gcContentArray[l]
                    next_gc = gcContentArray[l+1]
                    p4, = ax4.plot((position,nextPosition),(gc,next_gc), color='#824C71', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                secondToLastPos = positionArray[-2]
                lastPos = positionArray[-1]
                secondToLastPos = secondToLastPos + (windowSize/2)
                secondToLast_gc = gcContentArray[-2]
                last_gc = gcContentArray[-1]
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPos = lastPos + (lastWindowLength/2)
                else:
                    lastPos = lastPos + (windowSize/2)
                ax4.plot((secondToLastPos,lastPos),(secondToLast_gc,last_gc), color='#824C71', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                ax4.spines["left"].set_visible(True)
                ax4.spines["right"].set_visible(False)
                ax4.yaxis.set_label_position('left')
                ax4.yaxis.set_ticks_position('left')
                ax4.spines['left'].set_color('#824C71')
                ax4.tick_params(axis='y', colors='#824C71')
                ax4.set_ylabel('GC content (%)')
                ax4.yaxis.get_label().set_color(p4.get_color())
                ax4.yaxis.set_tick_params(labelright=False)
                ax4.margins(0)
                ax4.set_ylim(0,100)
                ticks = ax4.get_xticks()/1000.0
                ax4.set_xticks(ticks)
                ax4.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='#824C71', marker='o',markersize=6, label='GC content (%)', alpha=1.0, fillstyle='full', markeredgecolor='#824C71', markeredgewidth=0.0, linestyle='dashed')
                if 'GC content (%)' not in legendDict1:
                    legendDict1['GC content (%)'] = legendInfo
            ax5 = ax1.twinx()
            ax5.grid(False)
            if assemblyID in geneDensity and chromID in geneDensity[assemblyID]:
                positions,counts = zip(*geneDensity[assemblyID][chromID])
                minCount = min(counts)
                maxCount = max(counts)
                geneArray = np.array(geneDensity[assemblyID][chromID])
                for k in range(len(geneArray)-2):
                    position,geneDensityValue = geneArray[k]
                    nextPosition,nextGeneDensityValue = geneArray[k+1]
                    position = position + (windowSize/2)
                    nextPosition = nextPosition + (windowSize/2)
                    p5, = ax5.plot((position,nextPosition),(geneDensityValue,nextGeneDensityValue), color='blue', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                secondToLastPos,secondToLastGeneDensityValue = geneArray[-2]
                secondToLastPos = secondToLastPos + (windowSize/2)
                lastPos,lastGeneDensityValue = geneArray[-1]
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPos = lastPos + (lastWindowLength/2)
                else:
                    lastPos = lastPos + (windowSize/2)
                ax5.plot((secondToLastPos,lastPos),(secondToLastGeneDensityValue,lastGeneDensityValue), color='blue', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                ax5.yaxis.set_tick_params(labelleft=False)
                ax5.yaxis.set_tick_params(labelright=False)
                ax5.spines["right"].set_visible(False)
                ax5.xaxis.set_ticks_position('none')
                ax5.set_yticks([])
                ax5.margins(0)
                ax5.set_ylim(0,)
                ticks = ax5.get_xticks()/1000.0
                ax5.set_xticks(ticks)
                ax5.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='blue', marker='o',markersize=6, label='Transcript count (' + str(minCount) + '-' + str(maxCount) + ' per ' + str(windowSize) + ' bp bins)', alpha=1.0, fillstyle='full', markeredgecolor='blue', markeredgewidth=0.0, linestyle='dashed')
                if 'Transcript count' not in legendDict1:
                    legendDict1['Transcript count'] = legendInfo
            #ax7.axis("off")
            legendList1 = list(legendDict1.values())
            #ax7.legend(handles=legendList1, frameon=False, ncol=4, loc='upper center')
            ax1.legend(handles=legendList1, frameon=False, ncol=3, bbox_to_anchor=(0.5, -0.5), loc='center')
            #fig.supxlabel('Chromosome position (kb)')
            ax1.set_xlabel('Chromosome position (kb)', fontsize=16)
            ax1.set_title(assemblyID + " " + chromID + " " + str(seqLen) + " bp\nGenomic window " + str(windowStart) + "-" + str(windowStop) + " bp, " + str(windowSize) + ' bp bins')
            plt.savefig(assemblyID + "_" + chromID + "_" + featureType + "_distribution_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + "_condensed.png",dpi=600)
            plt.savefig(assemblyID + "_" + chromID + "_" + featureType + "_distribution_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + "_condensed.svg",)
            plt.close()

                    
def plotData(densityData,compiledTimes,methylationDensity,chromosomeLengthDict,gcDataDict,geneDensity,lastWindowRemainderDict,featureType,windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict,ty1_specificChromDataDict,ty3_specificChromDataDict,teColorMap,ltrColorMap,windowStart,windowStop,synthaseCoordDict,alt4CoordDict,bkrCoordDict,specificTE):
    windowSizeThreshold = 1000000
    for assemblyID in densityData:
        for chromID in densityData[assemblyID]:
            fullChromID = assemblyID + '.' + chromID
            seqLen = chromosomeLengthDict[assemblyID][chromID]
            lastWindowLength = lastWindowRemainderDict[assemblyID][chromID]
            if assemblyID in dup_variantDensityDict and chromID in dup_variantDensityDict[assemblyID] and assemblyID in ty1_specificChromDataDict and chromID in ty1_specificChromDataDict[assemblyID] and assemblyID in ty3_specificChromDataDict and chromID in ty3_specificChromDataDict[assemblyID]:
                fig, (ax1, ax7, ax8, ax10) = plt.subplots(nrows=4, sharex=True, figsize=(25,12.5), gridspec_kw=dict(height_ratios=[3, 1.5, 1, 1.5]))
            else:
                fig, ax1 = plt.subplots(figsize=(25,12.5))
            plt.rcParams['font.size'] = 16
            legendDict1 = {}
            legendDict2 = {}
            legendDict3 = {}
            legendDict4 = {}
            for featureID in densityData[assemblyID][chromID]:
                featureArray = np.array(densityData[assemblyID][chromID][featureID])
                if featureType == 'LTRs':
                    colorMap = ltrColorMap
                else:
                    colorMap = teColorMap
                for i in range(len(featureArray)-1):
                    position,teCount,teTotalLen,totalTEPercent = featureArray[i]
                    ax1.bar(position,totalTEPercent,windowSize, align='edge', label=featureID, alpha=1, edgecolor=colorMap[featureID], color='None', linewidth=4)
                lastPosition,lastTECount,lastTotalTELen,lastTotalTEPercent = featureArray[-1]
                # this is assessing whether the last window is less than the span of a full window length -- otherwise, it's the length that is leftover
                testStop = lastPosition+windowSize
                if testStop > seqLen:
                    ax1.bar(lastPosition,lastTotalTEPercent,lastWindowLength, align='edge', label=featureID, alpha=1, edgecolor=colorMap[featureID], color='None', linewidth=4)
                else:
                    ax1.bar(lastPosition,lastTotalTEPercent,windowSize, align='edge', label=featureID, alpha=1, edgecolor=colorMap[featureID], color='None', linewidth=4)
                ax1.spines['left'].set_color('black')
                ax1.spines['right'].set_color('black')
                ax1.tick_params(axis='y', colors='black')
                ax1.set_ylim(0,100)
                ax1.set_ylabel('Transposon percent (%)', fontsize=16)
                ax1.margins(0)
                legendInfo = mlines.Line2D([], [], color=colorMap[featureID], marker='s',linestyle="None",markersize=10, label=featureID, alpha=1.0, fillstyle='full', markeredgecolor=colorMap[featureID], markeredgewidth=0.0)
                ticks = ax1.get_xticks()/1000.0
                ax1.set_xticks(ticks)
                ax1.set_xticklabels(ticks.astype(int))
                if featureID not in legendDict1:
                    legendDict1[featureID] = legendInfo
                    
            synthaseColors = {'THCAS':'blue', 'CBDAS':'red', 'CBCAS':'orange'}
            if fullChromID in synthaseCoordDict:
                for moleculeType in synthaseCoordDict[fullChromID]:
                    for synthaseStart,synthaseStop,synthaseID in synthaseCoordDict[fullChromID][moleculeType]:
                        if synthaseStart >= windowStart and synthaseStart < windowStop:
                            if windowSize >= windowSizeThreshold:
                                ax1.vlines(synthaseStart, 0, 100, linewidth=4, color=synthaseColors[moleculeType])
                            else:
                                ax1.axvspan(synthaseStart,synthaseStop, ymin=0, ymax=100, color=synthaseColors[moleculeType], alpha=1, linestyle='solid', linewidth=2)
                            if len(synthaseCoordDict[fullChromID][moleculeType]) > 1:
                                synthaseLegendInfo = mlines.Line2D([], [], color=synthaseColors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + ', ' + str(len(synthaseCoordDict[fullChromID][moleculeType])) + ' copies', alpha=1.0, fillstyle='full', markeredgecolor=synthaseColors[moleculeType], markeredgewidth=0.0, linewidth=4)
                            else:
                                synthaseLegendInfo = mlines.Line2D([], [], color=synthaseColors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + '\n' + str(synthaseStart) + '-' + str(synthaseStop) + ' bp', alpha=1.0, fillstyle='full', markeredgecolor=synthaseColors[moleculeType], markeredgewidth=0.0, linewidth=4)
                            if moleculeType not in legendDict1:
                                legendDict1[moleculeType] = synthaseLegendInfo

            alt4Colors = {'ALT4':'purple'}
            if fullChromID in alt4CoordDict:
                for moleculeType in alt4CoordDict[fullChromID]:
                    for alt4Start,alt4Stop,alt4ID in alt4CoordDict[fullChromID][moleculeType]:
                        if alt4Start >= windowStart and alt4Start < windowStop:
                            if windowSize >= windowSizeThreshold:
                                ax1.vlines(alt4Start, 0, 100, linewidth=4, color=alt4Colors[moleculeType])
                            else:
                                ax1.axvspan(alt4Start,alt4Stop, ymin=0, ymax=100, color=alt4Colors[moleculeType], alpha=1, linestyle='solid', linewidth=2)
                            if len(alt4CoordDict[fullChromID][moleculeType]) > 1:
                                alt4LegendInfo = mlines.Line2D([], [], color=alt4Colors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + ', ' + str(len(alt4CoordDict[fullChromID][moleculeType])) + ' copies', alpha=1.0, fillstyle='full', markeredgecolor=alt4Colors[moleculeType], markeredgewidth=0.0, linewidth=4)
                            else:
                                alt4LegendInfo = mlines.Line2D([], [], color=alt4Colors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + '\n' + str(alt4Start) + '-' + str(alt4Stop) + ' bp', alpha=1.0, fillstyle='full', markeredgecolor=alt4Colors[moleculeType], markeredgewidth=0.0, linewidth=4)
                            if moleculeType not in legendDict1:
                                legendDict1[moleculeType] = alt4LegendInfo
            bkrColors = {'BKR':'yellow'}
            if fullChromID in bkrCoordDict:
                for moleculeType in bkrCoordDict[fullChromID]:
                    for bkrStart,bkrStop,bkrID in bkrCoordDict[fullChromID][moleculeType]:
                        if bkrStart >= windowStart and bkrStart < windowStop:
                            if windowSize >= windowSizeThreshold:
                                ax1.vlines(bkrStart, 0, 100, linewidth=4, color=bkrColors[moleculeType])
                            else:
                                ax1.axvspan(bkrStart,bkrStop, ymin=0, ymax=100, color=bkrColors[moleculeType], alpha=1, linestyle='solid', linewidth=2)
                            if len(bkrCoordDict[fullChromID][moleculeType]) > 1:
                                bkrLegendInfo = mlines.Line2D([], [], color=bkrColors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + ', ' + str(len(bkrCoordDict[fullChromID][moleculeType])) + ' copies', alpha=1.0, fillstyle='full', markeredgecolor=bkrColors[moleculeType], markeredgewidth=0.0, linewidth=4)
                            else:
                                bkrLegendInfo = mlines.Line2D([], [], color=bkrColors[moleculeType], marker='None',linestyle="solid",markersize=10, label=moleculeType + '\n' + str(bkrStart) + '-' + str(bkrStop) + ' bp', alpha=1.0, fillstyle='full', markeredgecolor=bkrColors[moleculeType], markeredgewidth=0.0, linewidth=4)
                            if moleculeType not in legendDict1:
                                legendDict1[moleculeType] = bkrLegendInfo
                                
            ax2 = ax1.twinx()
            if assemblyID in methylationDensity and chromID in methylationDensity[assemblyID]:
                methylationArray = np.array(methylationDensity[assemblyID][chromID])
                for j in range(len(methylationArray)-2):
                    position,methylationProb = methylationArray[j]
                    # this position modification is to plot the marker in the center of the window
                    position = position + (windowSize/2)
                    nextPosition,nextMethylationProb = methylationArray[j+1]
                    nextPosition = nextPosition + (windowSize/2)
                    ax2.plot((position,nextPosition),(methylationProb,nextMethylationProb), color='black', marker='o', linewidth=4, markersize=6)
                # define last position explicitly
                # this position modification is to plot the marker in the center of the window 
                secondToLastPosition,secondToLastMethylationProb = methylationArray[-2]
                lastPos,lastMethylationProb = methylationArray[-1]
                secondToLastPosition = secondToLastPosition + (windowSize/2)
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPosition = lastPos + (lastWindowLength/2)
                else:
                    lastPosition = lastPos + (windowSize/2)
                ax2.plot((secondToLastPosition,lastPosition),(secondToLastMethylationProb,lastMethylationProb), color='black', marker='o', linewidth=4, markersize=6)
                ax2.set_ylabel('Average methylation probability')
                ax2.margins(0)
                ax2.set_ylim(0,100)
                ticks = ax2.get_xticks()/1000.0
                ax2.set_xticks(ticks)
                ax2.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='black', marker='o',markersize=6, label='Methylation probability', alpha=1.0, fillstyle='full', markeredgecolor='black', markeredgewidth=0.0)
                if 'Methylation_probability' not in legendDict1:
                    legendDict1['Methylation_probability'] = legendInfo

            ax3 = ax1.twinx()
            ax3.grid(False)
            ax3.spines.right.set_position(("axes", 1.08))
            if assemblyID in gcDataDict and chromID in gcDataDict[assemblyID]:
                positions,observed_expected,gcContent = zip(*gcDataDict[assemblyID][chromID])
                smallestPos = min(positions)
                largestPos = max(positions)
                largestObservedExpected = max(observed_expected)
                positionArray = np.array(positions)
                obs_expArray = np.array(observed_expected)
                for k in range(len(obs_expArray)-2):
                    position = positionArray[k]
                    # this position modification is to plot the marker in the center of the window
                    position = position + (windowSize/2)
                    obs_exp = obs_expArray[k]
                    nextPosition = positionArray[k+1]
                    nextPosition = nextPosition + (windowSize/2)
                    next_obs_exp = obs_expArray[k+1]
                    p3, = ax3.plot((position,nextPosition),(obs_exp,next_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=6)
                # this position modification is to plot the marker in the center of the window
                secondToLastPos = positionArray[-2]
                lastPos = positionArray[-1]
                secondToLastPos = secondToLastPos + (windowSize/2)
                secondToLast_obs_exp = obs_expArray[-2]
                last_obs_exp = obs_expArray[-1]
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPos = lastPos + (lastWindowLength/2)
                    ax3.hlines(y=obs_exp_threshold, xmin=smallestPos, xmax=largestPos+lastWindowLength, linewidth=3, color='#E0A458')
                else:
                    lastPos = lastPos + (windowSize/2)
                    ax3.hlines(y=obs_exp_threshold, xmin=smallestPos, xmax=largestPos+windowSize, linewidth=3, color='#E0A458')
                ax3.plot((secondToLastPos,lastPos),(secondToLast_obs_exp,last_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=6)
                ax3.spines['right'].set_color('#E0A458')
                ax3.tick_params(axis='y', colors='#E0A458')
                ax3.set_ylabel('Observed/expected CpG')
                ax3.yaxis.get_label().set_color(p3.get_color())
                ax3.margins(0)
                if largestObservedExpected > 1:
                    ax3.set_ylim(0,)
                else:
                    ax3.set_ylim(0,1)
                ticks = ax3.get_xticks()/1000.0
                ax3.set_xticks(ticks)
                ax3.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='#E0A458', marker='o',markersize=6, label='Observed/expected CpG', alpha=1.0, fillstyle='full', markeredgecolor='#E0A458', markeredgewidth=0.0)
                if 'Observed/expected CpG' not in legendDict1:
                    legendDict1['Observed/expected CpG'] = legendInfo

            ax4 = ax1.twinx()
            ax4.grid(False)
            ax4.spines.left.set_position(("axes", -0.08))
            if assemblyID in gcDataDict and chromID in gcDataDict[assemblyID]:
                positions,observed_expected,gcContent = zip(*gcDataDict[assemblyID][chromID])
                positionArray = np.array(positions)
                gcContentArray = np.array(gcContent)
                for l in range(len(gcContentArray)-2):
                    position = positionArray[l]
                    position = position + (windowSize/2)
                    nextPosition = positionArray[l+1]
                    nextPosition = nextPosition + (windowSize/2)
                    gc = gcContentArray[l]
                    next_gc = gcContentArray[l+1]
                    p4, = ax4.plot((position,nextPosition),(gc,next_gc), color='#824C71', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                secondToLastPos = positionArray[-2]
                lastPos = positionArray[-1]
                secondToLastPos = secondToLastPos + (windowSize/2)
                secondToLast_gc = gcContentArray[-2]
                last_gc = gcContentArray[-1]
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPos = lastPos + (lastWindowLength/2)
                else:
                    lastPos = lastPos + (windowSize/2)
                ax4.plot((secondToLastPos,lastPos),(secondToLast_gc,last_gc), color='#824C71', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                ax4.spines["left"].set_visible(True)
                ax4.spines["right"].set_visible(False)
                ax4.yaxis.set_label_position('left')
                ax4.yaxis.set_ticks_position('left')
                ax4.spines['left'].set_color('#824C71')
                ax4.tick_params(axis='y', colors='#824C71')
                ax4.set_ylabel('GC content (%)')
                ax4.yaxis.get_label().set_color(p4.get_color())
                ax4.yaxis.set_tick_params(labelright=False)
                ax4.margins(0)
                ax4.set_ylim(0,100)
                ticks = ax4.get_xticks()/1000.0
                ax4.set_xticks(ticks)
                ax4.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='#824C71', marker='o',markersize=6, label='GC content (%)', alpha=1.0, fillstyle='full', markeredgecolor='#824C71', markeredgewidth=0.0, linestyle='dashed')
                if 'GC content (%)' not in legendDict1:
                    legendDict1['GC content (%)'] = legendInfo

            ax5 = ax1.twinx()
            ax5.grid(False)
            if assemblyID in geneDensity and chromID in geneDensity[assemblyID]:
                positions,counts = zip(*geneDensity[assemblyID][chromID])
                minCount = min(counts)
                maxCount = max(counts)
                geneArray = np.array(geneDensity[assemblyID][chromID])
                for k in range(len(geneArray)-2):
                    position,geneDensityValue = geneArray[k]
                    nextPosition,nextGeneDensityValue = geneArray[k+1]
                    position = position + (windowSize/2)
                    nextPosition = nextPosition + (windowSize/2)
                    p5, = ax5.plot((position,nextPosition),(geneDensityValue,nextGeneDensityValue), color='blue', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                secondToLastPos,secondToLastGeneDensityValue = geneArray[-2]
                secondToLastPos = secondToLastPos + (windowSize/2)
                lastPos,lastGeneDensityValue = geneArray[-1]
                testStop = lastPos+windowSize
                if testStop > seqLen:
                    lastPos = lastPos + (lastWindowLength/2)
                else:
                    lastPos = lastPos + (windowSize/2)
                ax5.plot((secondToLastPos,lastPos),(secondToLastGeneDensityValue,lastGeneDensityValue), color='blue', marker='o', linewidth=4, markersize=6, linestyle='dashed')
                ax5.yaxis.set_tick_params(labelleft=False)
                ax5.yaxis.set_tick_params(labelright=False)
                ax5.spines["right"].set_visible(False)
                ax5.xaxis.set_ticks_position('none')
                ax5.set_yticks([])
                ax5.margins(0)
                ax5.set_ylim(0,)
                ticks = ax5.get_xticks()/1000.0
                ax5.set_xticks(ticks)
                ax5.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='blue', marker='o',markersize=6, label='Transcript count (' + str(minCount) + '-' + str(maxCount) + ' per ' + str(windowSize) + ' bp bins)', alpha=1.0, fillstyle='full', markeredgecolor='blue', markeredgewidth=0.0, linestyle='dashed')
                if 'Transcript count' not in legendDict1:
                    legendDict1['Transcript count'] = legendInfo

            ax7.axis("off")
            legendList1 = list(legendDict1.values())
            ax7.legend(handles=legendList1, frameon=False, ncol=4, loc='upper center')
            variantLabelMap = {'dups':'Duplication', 'transloc':'Translocation', 'inv':'Inversion'}
            variantColorMap = {'dups':'purple', 'transloc':'red', 'inv':'orange'}

            if assemblyID in dup_variantDensityDict and chromID in dup_variantDensityDict[assemblyID]:
                dupArray = np.array(dup_variantDensityDict[assemblyID][chromID])
                for m in range(len(dupArray)-2):
                    dup_position,dup_variantCount,dup_totalVariantLength,dup_totalVariantPercent = dupArray[m]
                    next_dup_position,next_dup_variantCount,next_dup_totalVariantLength,next_dup_totalVariantPercent = dupArray[m+1]
                    dup_position = dup_position + (windowSize/2)
                    next_dup_position = next_dup_position + (windowSize/2)
                    dup_p8, = ax8.plot((dup_position,next_dup_position),(dup_totalVariantPercent,next_dup_totalVariantPercent), ls='-', color=variantColorMap['dups'], label=variantLabelMap['dups'], linewidth=4)
                dup_secondToLastPos,dup_secondToLastVariantCount,dup_secondToLastTotalVariantLength,dup_secondToLastTotalVariantPercent = dupArray[-2]
                dup_secondToLastPos = dup_secondToLastPos + (windowSize/2)
                dup_lastPos,dup_lastVariantCount,dup_lastTotalVariantLength,dup_lastTotalVariantPercent = dupArray[-1]
                testStop = dup_lastPos+windowSize
                if testStop > seqLen:
                    dup_lastPos = dup_lastPos + (lastWindowLength/2)
                else:
                    dup_lastPos = dup_lastPos + (windowSize/2)
                ax8.plot((dup_secondToLastPos,dup_lastPos),(dup_secondToLastTotalVariantPercent,dup_lastTotalVariantPercent), ls='-', color=variantColorMap['dups'], label=variantLabelMap['dups'], linewidth=4)
                dupLegendInfo1 = mlines.Line2D([], [], color=variantColorMap['dups'], marker='o', markersize=10, label=variantLabelMap['dups'] + ' percent', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['dups'], markeredgewidth=0.0)
                if 'Duplication percent' not in legendDict2:
                    legendDict2['Duplication percent'] = dupLegendInfo1

            if assemblyID in inv_variantDensityDict and chromID in inv_variantDensityDict[assemblyID]:
                invArray = np.array(inv_variantDensityDict[assemblyID][chromID])
                for m in range(len(invArray)-2):
                    inv_position,inv_variantCount,inv_totalVariantLength,inv_totalVariantPercent = invArray[m]
                    next_inv_position,next_inv_variantCount,next_inv_totalVariantLength,next_inv_totalVariantPercent = invArray[m+1]
                    inv_position = inv_position + (windowSize/2)
                    next_inv_position = next_inv_position + (windowSize/2)
                    inv_p8, = ax8.plot((inv_position,next_inv_position),(inv_totalVariantPercent,next_inv_totalVariantPercent), ls='-', color=variantColorMap['inv'], label=variantLabelMap['inv'], linewidth=4)
                inv_secondToLastPos,inv_secondToLastVariantCount,inv_secondToLastTotalVariantLength,inv_secondToLastTotalVariantPercent = invArray[-2]
                inv_secondToLastPos = inv_secondToLastPos + (windowSize/2)
                inv_lastPos,inv_lastVariantCount,inv_lastTotalVariantLength,inv_lastTotalVariantPercent = invArray[-1]
                testStop = inv_lastPos+windowSize
                if testStop > seqLen:
                    inv_lastPos = inv_lastPos + (lastWindowLength/2)
                else:
                    inv_lastPos = inv_lastPos + (windowSize/2)
                ax8.plot((inv_secondToLastPos,inv_lastPos),(inv_secondToLastTotalVariantPercent,inv_lastTotalVariantPercent), ls='-', color=variantColorMap['inv'], label=variantLabelMap['inv'], linewidth=4)
                invLegendInfo1 = mlines.Line2D([], [], color=variantColorMap['inv'], marker='o', markersize=10, label=variantLabelMap['inv'] + ' percent', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['inv'], markeredgewidth=0.0)
                if 'Inversion percent' not in legendDict2:
                    legendDict2['Inversion percent'] = invLegendInfo1
                
            if assemblyID in trans_variantDensityDict and chromID in trans_variantDensityDict[assemblyID]:
                transArray = np.array(trans_variantDensityDict[assemblyID][chromID])
                for m in range(len(transArray)-2):
                    trans_position,trans_variantCount,trans_totalVariantLength,trans_totalVariantPercent = transArray[m]
                    next_trans_position,next_trans_variantCount,next_trans_totalVariantLength,next_trans_totalVariantPercent = transArray[m+1]
                    trans_position = trans_position + (windowSize/2)
                    next_trans_position = next_trans_position + (windowSize/2)
                    trans_p8, = ax8.plot((trans_position,next_trans_position),(trans_totalVariantPercent,next_trans_totalVariantPercent), ls='-', color=variantColorMap['transloc'], label=variantLabelMap['transloc'], linewidth=4)
                trans_secondToLastPos,trans_secondToLastVariantCount,trans_secondToLastTotalVariantLength,trans_secondToLastTotalVariantPercent = transArray[-2]
                trans_secondToLastPos = trans_secondToLastPos + (windowSize/2)
                trans_lastPos,trans_lastVariantCount,trans_lastTotalVariantLength,trans_lastTotalVariantPercent = transArray[-1]
                testStop = trans_lastPos+windowSize
                if testStop > seqLen:
                    trans_lastPos = trans_lastPos + (lastWindowLength/2)
                else:
                    trans_lastPos = trans_lastPos + (windowSize/2)
                ax8.plot((trans_secondToLastPos,trans_lastPos),(trans_secondToLastTotalVariantPercent,trans_lastTotalVariantPercent), ls='-', color=variantColorMap['transloc'], label=variantLabelMap['transloc'], linewidth=4)
                transLegendInfo1 = mlines.Line2D([], [], color=variantColorMap['transloc'], marker='o', markersize=10, label=variantLabelMap['transloc'] + ' percent', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['transloc'], markeredgewidth=0.0)
                if 'Translocation percent' not in legendDict2:
                    legendDict2['Translocation percent'] = transLegendInfo1

            ax8.set_ylabel('Variant percent (%)', rotation=45, horizontalalignment='right', fontsize=16)
            ax8.set_ylim(0,100)
            ax8.margins(0)
            ticks = ax8.get_xticks()/1000.0
            ax8.set_xticks(ticks)
            ax8.set_xticklabels(ticks.astype(int))

            if 'Ty1' in ty1_specificChromDataDict[assemblyID][chromID]:
                ltrTypeLabel = 'Ty1_LTR_retrotransposon'
                ty1Array = np.array(ty1_specificChromDataDict[assemblyID][chromID]['Ty1'])
                positions,soloCounts,intactCounts,siRatios = zip(*ty1_specificChromDataDict[assemblyID][chromID]['Ty1'])
                maxSoloCount = max(soloCounts)
                maxYValue = max(max(soloCounts),max(intactCounts))
                if maxYValue == 0:
                    maxYValue = 1.0
                ty1Array = np.array(ty1_specificChromDataDict[assemblyID][chromID]['Ty1'])
                for i in range(len(ty1Array)-1):
                    position,soloCount,intactCount,siRatio = ty1Array[i]
                    position = int(position)
                    siRatio = float(siRatio)
                    ax10.bar(position,siRatio,windowSize, align='edge', label='Solo:intact ratio', edgecolor=ltrColorMap[ltrTypeLabel], color='None', linewidth=4)
                last_position,last_soloCount,last_intactCount,last_siRatio = ty1Array[-1]
                last_position = int(last_position)
                last_siRatio = float(last_siRatio)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    ax10.bar(last_position,last_siRatio,lastWindowLength, align='edge', edgecolor=ltrColorMap[ltrTypeLabel], color='None', linewidth=4)
                else:
                    ax10.bar(last_position,last_siRatio,windowSize, align='edge', edgecolor=ltrColorMap[ltrTypeLabel], color='None', linewidth=4)
                ax10.set_ylabel('Average solo:intact\nLTR ratio', rotation=45, horizontalalignment='right', fontsize=16)
                ax10.margins(0)
                ax10.set_ylim(0,)
                ticks = ax10.get_xticks()/1000.0
                ax10.set_xticks(ticks)
                ax10.set_xticklabels(ticks.astype(int))

                ax10.xaxis.set_tick_params(labelsize=16)
                
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='s',linestyle="None",markersize=10, label='Ty1-LTR solo:intact ratio', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                if 'Ty1 solo:intact ratio' not in legendDict4:
                    legendDict4['Ty1 solo:intact ratio'] = legendInfo
            # plot solo and intact Ty3-LTR info
            if 'Ty3' in ty3_specificChromDataDict[assemblyID][chromID]:
                #ax13 = ax10.twinx()
                ltrTypeLabel = 'Ty3_LTR_retrotransposon'
                ty3Array = np.array(ty3_specificChromDataDict[assemblyID][chromID]['Ty3'])
                positions,soloCounts,intactCounts,siRatios = zip(*ty3_specificChromDataDict[assemblyID][chromID]['Ty3'])
                maxSoloCount = max(soloCounts)
                maxYValue = max(max(soloCounts),max(intactCounts))
                if maxYValue == 0:
                    maxYValue = 1.0
                ty3Array = np.array(ty3_specificChromDataDict[assemblyID][chromID]['Ty3'])
                for i in range(len(ty3Array)-1):
                    position,soloCount,intactCount,siRatio = ty3Array[i]
                    position = int(position)
                    siRatio = float(siRatio)
                    ax10.bar(position,siRatio,windowSize, align='edge', label='Solo:intact ratio', edgecolor=ltrColorMap[ltrTypeLabel], color='None', linewidth=4)
                last_position,last_soloCount,last_intactCount,last_siRatio = ty3Array[-1]
                last_position = int(last_position)
                last_siRatio = float(last_siRatio)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    ax10.bar(last_position,last_siRatio,lastWindowLength, align='edge', edgecolor=ltrColorMap[ltrTypeLabel], color='None', linewidth=4)
                else:
                    ax10.bar(last_position,last_siRatio,windowSize, align='edge', edgecolor=ltrColorMap[ltrTypeLabel], color='None', linewidth=4)
                ax10.margins(0)
                ax10.set_ylim(0,)
                ticks = ax10.get_xticks()/1000.0
                ax10.set_xticks(ticks)
                ax10.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='s',linestyle="None",markersize=10, label='Ty3-LTR solo:intact ratio', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                if 'Ty3 solo:intact ratio' not in legendDict4:
                    legendDict4['Ty3 solo:intact ratio'] = legendInfo
            legendList4 = list(legendDict4.values())
            ax10.legend(handles=legendList4, frameon=False, ncol=1, loc='upper right')
            ax10.set_xlabel('Chromosome position (kb)', fontsize=20)
            legendList2 = list(legendDict2.values())
            legendList3 = list(legendDict3.values())
            ax8.legend(handles=legendList2, frameon=False, ncol=1, loc='upper right') 
            ax1.set_title(assemblyID + " " + chromID + " " + str(seqLen) + " bp\nGenomic window " + str(windowStart) + "-" + str(windowStop) + " bp, " + str(windowSize) + ' bp bins')
            plt.savefig(assemblyID + "_" + chromID + "_" + featureType + "_distribution_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".png",dpi=600)
            plt.savefig(assemblyID + "_" + chromID + "_" + featureType + "_distribution_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".svg",)
            plt.close()


def plotSimplifiedTimeData(compiledTimes,chromosomeLengthDict,colorMap,lastWindowRemainderDict,featureType,windowSize,windowStart,windowStop):
    colorDict = {'0_5':'#AA3377', '6_10':'#EE6677', '11_15':'#66CCEE', '16_20':'#228833', '21_25':'#4477AA'}
    labelDict = {'0_5':'<6 mya', '6_10':'>=6 to <11 mya', '11_15':'>=11 to <16 mya', '16_20':'>=16 to <21 mya', '21_25':'>=21 to <26 mya'}
    for assemblyID in compiledTimes:
        for chromID in compiledTimes[assemblyID]:
            seqLen = chromosomeLengthDict[assemblyID][chromID]
            lastWindowLength = lastWindowRemainderDict[assemblyID][chromID]
            for featureID in compiledTimes[assemblyID][chromID]:
                fig, ax1 = plt.subplots(figsize=(20,5))
                plt.rcParams['font.size'] = 16
                legendDict = {}
                for timeID in compiledTimes[assemblyID][chromID][featureID]:
                    featureArray = np.array(compiledTimes[assemblyID][chromID][featureID][timeID])
                    for i in range(len(featureArray)-1):
                        position,density,totalTELen,totalTEPercent,avgTime,specificTimePercent = featureArray[i]
                        ax1.bar(position,specificTimePercent,windowSize, align='edge', label=featureID, alpha=0.75, color=colorDict[timeID])
                    ax1.spines['left'].set_color('black')
                    ax1.spines['right'].set_color('black')
                    ax1.tick_params(axis='y', colors='black')
                    plt.xlabel('Chromosome position (kb)')
                    plt.ylabel('Transposon percent (%)')
                    ax1.margins(0)
                    ax1.set_ylim(0,100)
                    legendInfo = mlines.Line2D([], [], color=colorDict[timeID], marker='s',linestyle="None",markersize=10, label=labelDict[timeID], alpha=1.0, fillstyle='full', markeredgecolor=colorDict[timeID], markeredgewidth=0.0)
                    ticks = ax1.get_xticks()/1000.0
                    ax1.set_xticks(ticks)
                    ax1.set_xticklabels(ticks.astype(int))
                    if timeID not in legendDict:
                        legendDict[timeID] = legendInfo

                legendList = list(legendDict.values())
                plt.title(assemblyID + " " + chromID + " " + featureID + " " + str(windowStart) + "-" + str(windowStop) + " bp")
                plt.legend(handles=legendList, frameon=False, ncol=2, loc='upper center')
                plt.tight_layout()
                plt.savefig(assemblyID + "_" + chromID + "_" + featureID + "_simplified_subset_time_densities_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".png",dpi=600)
                plt.savefig(assemblyID + "_" + chromID + "_" + featureID + "_simplified_subset_time_densities_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".svg")
                plt.close()


                
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
                
            
usage = "Usage: " + sys.argv[0] + " <gene data> <DNA transposon density/count data> <DNA transposon time data> <long terminal retrotransposon density/count data> <long terminal retrotransposon time data> <methylation data> <GC data> <duplication variant data> <inversion variant data> <translocation variant data> <chromosome lengths> <Ty1 solo/intact counts file> <Ty3 solo/intact counts file> <synthase bed file> <alt4 annotations> <bkr annotations> <orientation file> <assembly ID> <window size, e.g. 1000000> <window start> <window stop> <specific TE category (TE categories can include one of the following entries): PIF_Harbinger_TIR_transposon, Mutator_TIR_transposon, CACTA_TIR_transposon, hAT_TIR_transposon, helitron, Tc1_Mariner_TIR_transposon> <specific LTR category (LTR categories can include one of the following entries): unknown_LTR_retrotransposon, Ty1_LTR_retrotransposon, Ty3_LTR_retrotransposon>"
if len(sys.argv) != 24:
    print(usage)
    sys.exit()

geneData = sys.argv[1]
teData = sys.argv[2]
teTimeData = sys.argv[3]
ltrData = sys.argv[4]
ltrTimeData = sys.argv[5]
methylationData = sys.argv[6]
gcData = sys.argv[7]
dupData = sys.argv[8]
invData = sys.argv[9]
translocData = sys.argv[10]
chromosomeLengths = sys.argv[11]
ty1_soloIntactRatioFile = sys.argv[12]
ty3_soloIntactRatioFile = sys.argv[13]
synthaseBedFile = sys.argv[14]
alt4BedFile = sys.argv[15]
bkrBedFile = sys.argv[16]
orientationFile = sys.argv[17]
assemblyID = sys.argv[18]
windowSize = sys.argv[19]
windowStart = sys.argv[20]
windowStop = sys.argv[21]
specificTE = sys.argv[22]
specificLTR = sys.argv[23]

windowSize = int(windowSize)
windowStart = int(windowStart)
windowStop = int(windowStop)
obs_exp_threshold = 0.6

# orientation file for annotation coords
orientationDict = readOrientationFile(orientationFile)

# chromosome lengths
chromosomeLengthDict = readChromosomeLengthData(chromosomeLengths)
lastWindowRemainderDict = calculateLastWindowLength(chromosomeLengthDict,windowSize)

# ALT and BKR annotations
alt4CoordDict = readOtherAnnotationBedFiles(alt4BedFile,'ALT4',orientationDict,chromosomeLengthDict)
bkrCoordDict = readOtherAnnotationBedFiles(bkrBedFile,'BKR',orientationDict,chromosomeLengthDict)

# synthase bed file
synthaseCoordDict = readBedFile(synthaseBedFile,orientationDict,chromosomeLengthDict)

# solo/intact counts
ty1_specificChromDataDict = readSoloIntactRatioFile(ty1_soloIntactRatioFile,windowStart,windowStop)
ty3_specificChromDataDict = readSoloIntactRatioFile(ty3_soloIntactRatioFile,windowStart,windowStop)

# gene counts
geneDensity = readGeneData(geneData,windowStart,windowStop)

# transposable elements densities/counts
teDensityData,teFeatureIDs = readTEDensityData(teData,windowStart,windowStop)
ltrDensityData,ltrFeatureIDs = readTEDensityData(ltrData,windowStart,windowStop)

# transposable element times
teTimeDataDict = readTETimeData(teTimeData,windowStart,windowStop)
ltrTimeDataDict = readTETimeData(ltrTimeData,windowStart,windowStop)

# gc data
gcDataDict = readGCData(gcData,windowStart,windowStop)

# methylation data
methylationDensity = readMethylationData(methylationData,windowStart,windowStop)

# structural variants
dup_variantDensityDict = readVariantData(dupData,windowStart,windowStop)
inv_variantDensityDict = readVariantData(invData,windowStart,windowStop)
trans_variantDensityDict = readVariantData(translocData,windowStart,windowStop)

# color maps for plotting
teColorMap = createColorMapDict(teFeatureIDs,'TEs')
ltrColorMap = createColorMapDict(ltrFeatureIDs,'LTRs')

# plot data
plotData(teDensityData,teTimeDataDict,methylationDensity,chromosomeLengthDict,gcDataDict,geneDensity,lastWindowRemainderDict,'TEs',windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict,ty1_specificChromDataDict,ty3_specificChromDataDict,teColorMap,ltrColorMap,windowStart,windowStop,synthaseCoordDict,alt4CoordDict,bkrCoordDict,specificTE)

plotData(ltrDensityData,ltrTimeDataDict,methylationDensity,chromosomeLengthDict,gcDataDict,geneDensity,lastWindowRemainderDict,'LTRs',windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict,ty1_specificChromDataDict,ty3_specificChromDataDict,teColorMap,ltrColorMap,windowStart,windowStop,synthaseCoordDict,alt4CoordDict,bkrCoordDict,specificLTR)

plotSimplifiedTimeData(teTimeDataDict,chromosomeLengthDict,teColorMap,lastWindowRemainderDict,'TEs',windowSize,windowStart,windowStop)
plotSimplifiedTimeData(ltrTimeDataDict,chromosomeLengthDict,ltrColorMap,lastWindowRemainderDict,'LTRs',windowSize,windowStart,windowStop)

plotCondensed(teDensityData,teTimeDataDict,methylationDensity,chromosomeLengthDict,gcDataDict,geneDensity,lastWindowRemainderDict,'TEs',windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict,ty1_specificChromDataDict,ty3_specificChromDataDict,teColorMap,ltrColorMap,windowStart,windowStop,synthaseCoordDict,specificTE)
