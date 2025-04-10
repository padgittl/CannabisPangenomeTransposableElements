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

def readSoloIntactRatioFile(soloIntactRatioFile,windowStart,windowStop,allData):
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
                if 'soloCount' not in allData:
                    allData['soloCount'] = []
                allData['soloCount'].append(soloCount)
                if 'intactCount' not in allData:
                    allData['intactCount'] = []
                allData['intactCount'].append(intactCount)
                if genomicPosition >= windowStart and genomicPosition < windowStop:
                    if siRatio != 'nan':
                        siRatio = float(siRatio)
                        if 'siRatio' not in allData:
                            allData['siRatio'] = []
                        allData['siRatio'].append(siRatio)
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
    return(specificChromDataDict,allData)


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


'''
allData['soloCount'].append(soloCount)
allData['intactCount'].appendOA(intactCount)
allData['siRatio'].append(siRatio)
globalSoloIntactRatioYMax,globalCountYMax

'''
def plotData(allData,chromosomeLengthDict,lastWindowRemainderDict,windowSize,ty1_specificChromDataDict,ty3_specificChromDataDict,ltrColorMap,windowStart,windowStop,globalSoloIntactRatioYMax,globalCountYMax,specificTE):
    # maxYValue = max(allData['soloCount'],allData['intactCount'])
    # maxIntact = max(allData['intactCount'])
    dataMaxYValue = max(max(allData['soloCount']),max(allData['intactCount']))
    maxYValue = globalCountYMax
    maxRatio = globalSoloIntactRatioYMax
    if dataMaxYValue > globalCountYMax:
        # print(maxYValue)
        maxYValue = dataMaxYValue
        print(dataMaxYValue,'larger than globalCountYMax')
        sys.exit()
    dataMaxRatio = max(allData['siRatio'])
    if dataMaxRatio > globalSoloIntactRatioYMax:
        maxRatio = dataMaxRatio
        print(dataMaxRatio,'larger than globalSoloIntactRatioYMax')
        sys.exit()
    for assemblyID in ty1_specificChromDataDict:
        for chromID in ty1_specificChromDataDict[assemblyID]:
            fullChromID = assemblyID + '.' + chromID
            seqLen = chromosomeLengthDict[assemblyID][chromID]
            lastWindowLength = lastWindowRemainderDict[assemblyID][chromID]
            # print("lastWindowLength,seqLen,fullChromID,lastWindowLength+windowSize")
            # print(lastWindowLength,seqLen,fullChromID,lastWindowLength+windowSize)
            # fig, (ax1, ax6, ax7, ax8, ax9, ax10, ax13)
            if assemblyID in ty3_specificChromDataDict and chromID in ty3_specificChromDataDict[assemblyID]:
                fig, (ax10, ax13) = plt.subplots(nrows=2, sharex=True, figsize=(20,10), gridspec_kw=dict(height_ratios=[1, 1]))
            else:
                fig, ax1 = plt.subplots(figsize=(20,10))
            plt.rcParams['font.size'] = 16
            legendDict4 = {}
            legendDict5 = {}
            colorMap = ltrColorMap

            # plot solo and intact Ty1-LTR info
            if 'Ty1' in ty1_specificChromDataDict[assemblyID][chromID]:
                ltrTypeLabel = 'Ty1_LTR_retrotransposon'
                ty1Array = np.array(ty1_specificChromDataDict[assemblyID][chromID]['Ty1'])
                positions,soloCounts,intactCounts,siRatios = zip(*ty1_specificChromDataDict[assemblyID][chromID]['Ty1'])
                '''
                maxSoloCount = max(soloCounts)
                maxYValue = max(max(soloCounts),max(intactCounts))
                if maxYValue == 0:
                    maxYValue = 1.0
                '''
                ty1Array = np.array(ty1_specificChromDataDict[assemblyID][chromID]['Ty1'])
                for i in range(len(ty1Array)-1):
                    position,soloCount,intactCount,siRatio = ty1Array[i]
                    position = int(position)
                    siRatio = float(siRatio)
                    ax10.bar(position,siRatio,windowSize, align='edge', label='Solo:intact ratio', alpha=0.75, color=ltrColorMap[ltrTypeLabel])
                last_position,last_soloCount,last_intactCount,last_siRatio = ty1Array[-1]
                last_position = int(last_position)
                last_siRatio = float(last_siRatio)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    ax10.bar(last_position,last_siRatio,lastWindowLength, align='edge', alpha=0.75, color=ltrColorMap[ltrTypeLabel])
                else:
                    ax10.bar(last_position,last_siRatio,windowSize, align='edge', alpha=0.75, color=ltrColorMap[ltrTypeLabel])
                # plt.xlabel('Chromosome position (Mb)')
                ax10.set_ylabel('Ty1\nsolo:intact ratio', rotation=45, horizontalalignment='right')
                ax10.margins(0)
                # ax10.set_ylim(0,)
                ax10.set_ylim(0,maxRatio)
                ticks = ax10.get_xticks()/1000.0
                ax10.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='s',linestyle="None",markersize=14, label='Solo:intact ratio', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                if 'Solo:intact ratio' not in legendDict4:
                    legendDict4['Solo:intact ratio'] = legendInfo

                # solo count
                ax11 = ax10.twinx()
                for i in range(len(ty1Array)-2):
                    position,soloCount,intactCount,siRatio = ty1Array[i]
                    position = int(position)
                    position = position + (windowSize/2)
                    soloCount = int(soloCount)
                    next_position,next_soloCount,next_intactCount,next_siRatio = ty1Array[i+1]
                    next_position = int(next_position)
                    next_position = next_position + (windowSize/2)
                    next_soloCount = int(next_soloCount)
                    ax11.plot((position,next_position),(soloCount,next_soloCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dotted')
                second_to_last_position,second_to_last_soloCount,second_to_last_intactCount,second_to_last_siRatio = ty1Array[-2]
                last_position,last_soloCount,last_intactCount,last_siRatio = ty1Array[-1]
                second_to_last_position = int(second_to_last_position)
                last_position = int(last_position)
                second_to_last_position = second_to_last_position + (windowSize/2)
                #last_position = last_position + (lastWindowLength/2)
                second_to_last_soloCount = float(second_to_last_soloCount)
                last_soloCount = float(last_soloCount)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    last_position = last_position + (lastWindowLength/2)
                else:
                    last_position = last_position + (windowSize/2)
                ax11.plot((second_to_last_position,last_position),(second_to_last_soloCount,last_soloCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dotted')
                ax11.margins(0)
                ax11.set_ylim(0,maxYValue)
                ticks = ax11.get_xticks()/1000.0
                ax11.set_xticklabels(ticks.astype(int))
                #ax11.yaxis.set_tick_params(labelleft=False)
                #ax11.yaxis.set_tick_params(labelright=False)
                #ax11.xaxis.set_ticks_position('none')
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='o',linestyle="dotted",markersize=14, label='Solo LTR count', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                if 'Solo LTR count' not in legendDict4:
                    legendDict4['Solo LTR count'] = legendInfo
                    
                ax12 = ax10.twinx()
                for i in range(len(ty1Array)-2):
                    position,soloCount,intactCount,siRatio = ty1Array[i]
                    position = int(position)
                    position = position + (windowSize/2)
                    intactCount = int(intactCount)
                    next_position,next_soloCount,next_intactCount,next_siRatio = ty1Array[i+1]
                    next_position = int(next_position)
                    next_position = next_position + (windowSize/2)
                    next_intactCount = int(next_intactCount)
                    ax12.plot((position,next_position),(intactCount,next_intactCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dashdot')
                second_to_last_position,second_to_last_soloCount,second_to_last_intactCount,second_to_last_siRatio = ty1Array[-2]
                last_position,last_soloCount,last_intactCount,last_siRatio = ty1Array[-1]
                second_to_last_position = int(second_to_last_position)
                last_position = int(last_position)
                second_to_last_position = second_to_last_position + (windowSize/2)
                #last_position = last_position + (lastWindowLength/2)
                second_to_last_intactCount = float(second_to_last_intactCount)
                last_intactCount = float(last_intactCount)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    last_position = last_position + (lastWindowLength/2)
                else:
                    last_position = last_position + (windowSize/2)
                ax12.plot((second_to_last_position,last_position),(second_to_last_intactCount,last_intactCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dotted')
                ax12.margins(0)
                ax12.set_ylim(0,maxYValue)
                ticks = ax12.get_xticks()/1000.0
                ax12.set_xticklabels(ticks.astype(int))
                ax12.yaxis.set_tick_params(labelleft=False)
                ax12.yaxis.set_tick_params(labelright=False)
                ax12.xaxis.set_ticks_position('none')
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='o',linestyle="dashdot",markersize=14, label='Intact LTR count', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                if 'Intact LTR count' not in legendDict4:
                    legendDict4['Intact LTR count'] = legendInfo
                legendList4 = list(legendDict4.values())
                ax10.legend(handles=legendList4, frameon=False, ncol=3, loc='upper right')

            # plot solo and intact Ty3-LTR info
            if 'Ty3' in ty3_specificChromDataDict[assemblyID][chromID]:
                ltrTypeLabel = 'Ty3_LTR_retrotransposon'
                ty3Array = np.array(ty3_specificChromDataDict[assemblyID][chromID]['Ty3'])
                positions,soloCounts,intactCounts,siRatios = zip(*ty3_specificChromDataDict[assemblyID][chromID]['Ty3'])
                '''
                maxSoloCount = max(soloCounts)
                maxYValue = max(max(soloCounts),max(intactCounts))
                if maxYValue ==	0:
                    maxYValue =	1.0
                '''
                ty3Array = np.array(ty3_specificChromDataDict[assemblyID][chromID]['Ty3'])
                for i in range(len(ty3Array)-1):
                    position,soloCount,intactCount,siRatio = ty3Array[i]
                    position = int(position)
                    siRatio = float(siRatio)
                    ax13.bar(position,siRatio,windowSize, align='edge', label='Solo:intact ratio', alpha=0.75, color=ltrColorMap[ltrTypeLabel])
                last_position,last_soloCount,last_intactCount,last_siRatio = ty3Array[-1]
                last_position = int(last_position)
                last_siRatio = float(last_siRatio)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    ax13.bar(last_position,last_siRatio,lastWindowLength, align='edge', alpha=0.75, color=ltrColorMap[ltrTypeLabel])
                else:
                    ax13.bar(last_position,last_siRatio,windowSize, align='edge', alpha=0.75, color=ltrColorMap[ltrTypeLabel])
                # plt.xlabel('Chromosome position (Mb)')
                ax13.set_ylabel('Ty3\nsolo:intact ratio', rotation=45, horizontalalignment='right')
                ax13.margins(0)
                #ax13.set_ylim(0,)
                ax13.set_ylim(0,maxRatio)
                ticks = ax13.get_xticks()/1000.0
                ax13.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='s',linestyle="None",markersize=14, label='Solo:intact ratio', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                if 'Solo:intact ratio' not in legendDict5:
                    legendDict5['Solo:intact ratio'] = legendInfo

                ax14 = ax13.twinx()
                for i in range(len(ty3Array)-2):
                    position,soloCount,intactCount,siRatio = ty3Array[i]
                    position = int(position)
                    position = position + (windowSize/2)
                    soloCount = int(soloCount)
                    next_position,next_soloCount,next_intactCount,next_siRatio = ty3Array[i+1]
                    next_position = int(next_position)
                    next_position = next_position + (windowSize/2)
                    next_soloCount = int(next_soloCount)
                    ax14.plot((position,next_position),(soloCount,next_soloCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dotted')
                second_to_last_position,second_to_last_soloCount,second_to_last_intactCount,second_to_last_siRatio = ty3Array[-2]
                last_position,last_soloCount,last_intactCount,last_siRatio = ty3Array[-1]
                second_to_last_position = int(second_to_last_position)
                last_position = int(last_position)
                second_to_last_position = second_to_last_position + (windowSize/2)
                #last_position = last_position + (lastWindowLength/2)
                second_to_last_soloCount = float(second_to_last_soloCount)
                last_soloCount = float(last_soloCount)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    last_position = last_position + (lastWindowLength/2)
                else:
                    last_position = last_position + (windowSize/2)
                ax14.plot((second_to_last_position,last_position),(second_to_last_soloCount,last_soloCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dotted')
                ax14.margins(0)
                ax14.set_ylim(0,maxYValue)
                ticks = ax14.get_xticks()/1000.0
                ax14.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='o',linestyle="dotted",markersize=14, label='Solo LTR count', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                ax14.set_ylabel('Solo and intact\nLTR counts', rotation=90, horizontalalignment='left')
                if 'Solo LTR count' not in legendDict5:
                    legendDict5['Solo LTR count'] = legendInfo
                    
                ax15 = ax13.twinx()
                for i in range(len(ty3Array)-2):
                    position,soloCount,intactCount,siRatio = ty3Array[i]
                    position = int(position)
                    position = position + (windowSize/2)
                    intactCount = int(intactCount)
                    next_position,next_soloCount,next_intactCount,next_siRatio = ty3Array[i+1]
                    next_position = int(next_position)
                    next_position = next_position + (windowSize/2)
                    next_intactCount = int(next_intactCount)
                    ax15.plot((position,next_position),(intactCount,next_intactCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dashdot')
                    
                second_to_last_position,second_to_last_soloCount,second_to_last_intactCount,second_to_last_siRatio = ty3Array[-2]
                last_position,last_soloCount,last_intactCount,last_siRatio = ty3Array[-1]
                second_to_last_position = int(second_to_last_position)
                last_position = int(last_position)
                second_to_last_position = second_to_last_position + (windowSize/2)
                #last_position = last_position + (lastWindowLength/2)
                second_to_last_intactCount = float(second_to_last_intactCount)
                last_intactCount = float(last_intactCount)
                testStop = last_position+windowSize
                if testStop > seqLen:
                    last_position = last_position + (lastWindowLength/2)
                else:
                    last_position = last_position + (windowSize/2)
                ax15.plot((second_to_last_position,last_position),(second_to_last_intactCount,last_intactCount), color=ltrColorMap[ltrTypeLabel], marker='o', linewidth=6, markersize=8, linestyle='dashdot')
                ax15.margins(0)
                ax15.set_ylim(0,maxYValue)
                ticks = ax15.get_xticks()/1000.0
                ax15.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color=ltrColorMap[ltrTypeLabel], marker='o',linestyle="dashdot",markersize=14, label='Intact LTR count', alpha=1.0, fillstyle='full', markeredgecolor=ltrColorMap[ltrTypeLabel], markeredgewidth=0.0)
                if 'Intact LTR count' not in legendDict5:
                    legendDict5['Intact LTR count'] = legendInfo
                legendList5 = list(legendDict5.values())
                ax13.legend(handles=legendList5, frameon=False, ncol=3, loc='upper right')
            
            ax10.set_title(assemblyID + " " + chromID + " " + str(seqLen) + " bp\nGenomic window " + str(windowStart) + "-" + str(windowStop) + " bp, " + str(windowSize) + ' bp bins')
            plt.savefig(assemblyID + "_" + chromID + "_soloIntactData_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".png",dpi=600)
            plt.savefig(assemblyID + "_" + chromID + "_soloIntactData_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".svg",)
            print(assemblyID + "_" + chromID + "_soloIntactData_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".png")
            print(assemblyID + "_" + chromID + "_soloIntactData_" + str(windowStart) + "_" + str(windowStop) + "_window" + str(windowSize) + ".svg")
            plt.close()

                
usage = "Usage: " + sys.argv[0] + " <chromosome lengths file> <Ty1 solo:intact> <Ty3 solo:intact> <assembly ID> <window size> <window start> <window stop> \n "
if len(sys.argv) != 8:
    print(usage)
    sys.exit()

chromosomeLengths = sys.argv[1]
ty1_soloIntactRatioFile = sys.argv[2]
ty3_soloIntactRatioFile = sys.argv[3]
assemblyID = sys.argv[4]
windowSize = sys.argv[5]
windowStart = sys.argv[6]
windowStop = sys.argv[7]

windowSize = int(windowSize)
windowStart = int(windowStart)
windowStop = int(windowStop)

allData = {}

# chromosome lengths
chromosomeLengthDict = readChromosomeLengthData(chromosomeLengths)
lastWindowRemainderDict = calculateLastWindowLength(chromosomeLengthDict,windowSize)

# solo/intact counts
ty1_specificChromDataDict,allData = readSoloIntactRatioFile(ty1_soloIntactRatioFile,windowStart,windowStop,allData)
ty3_specificChromDataDict,allData = readSoloIntactRatioFile(ty3_soloIntactRatioFile,windowStart,windowStop,allData)

# color maps for plotting
ltrColorMap = {'Ty1_LTR_retrotransposon':'#542788', 'unknown_LTR_retrotransposon':'#998ec3', 'Ty3_LTR_retrotransposon':'#b35806'}

# plot data
globalSoloIntactRatioYMax = 40
globalCountYMax = 40
plotData(allData,chromosomeLengthDict,lastWindowRemainderDict,windowSize,ty1_specificChromDataDict,ty3_specificChromDataDict,ltrColorMap,windowStart,windowStop,globalSoloIntactRatioYMax,globalCountYMax,'LTRs')


