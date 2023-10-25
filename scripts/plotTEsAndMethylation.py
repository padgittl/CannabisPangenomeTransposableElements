import sys, re, os, math, statistics
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


########
# MAIN #
########


def readChromosomeLengthData(chromosomeLengths):
    chromosomeLengthDict = {}
    with open(chromosomeLengths,'r') as F:
        for line in F:
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


def readGeneData(geneData):
    geneDensity = {}
    with open(geneData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,position,geneCount = line.strip().split('\t')
                if assemblyID not in geneDensity:
                    geneDensity[assemblyID] = {}
                if chromID not in geneDensity[assemblyID]:
                    geneDensity[assemblyID][chromID] = []
                geneDensity[assemblyID][chromID].append((float(position),int(geneCount)))
    return(geneDensity)


def readTETimeData(timeData):
    timeDataDict = {}
    with open(timeData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,featureID,timeID,chromStartPosition,teCount,totalTELenPerWindow,percentWindow,averageTime,specificTimeTEPercent = line.strip().split('\t')
                if assemblyID not in timeDataDict:
                    timeDataDict[assemblyID] = {}
                if chromID not in timeDataDict[assemblyID]:
                    timeDataDict[assemblyID][chromID] = {}
                if featureID not in timeDataDict[assemblyID][chromID]:
                    timeDataDict[assemblyID][chromID][featureID] = {}
                if timeID not in timeDataDict[assemblyID][chromID][featureID]:
                    timeDataDict[assemblyID][chromID][featureID][timeID] = []
                timeDataDict[assemblyID][chromID][featureID][timeID].append((float(chromStartPosition),int(teCount),int(totalTELenPerWindow),float(percentWindow),float(averageTime),float(specificTimeTEPercent)))
    return(timeDataDict)
                

def readTEDensityData(teData):
    densityData = {}
    featureIDs = {}
    with open(teData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,featureID,position,teCount,totalTELenPerWindow,percentWindow = line.strip().split('\t')
                if featureID not in featureIDs:
                    featureIDs[featureID] = 1
                if assemblyID not in densityData:
                    densityData[assemblyID] = {}
                if chromID not in densityData[assemblyID]:
                    densityData[assemblyID][chromID] = {}
                if featureID not in densityData[assemblyID][chromID]:
                    densityData[assemblyID][chromID][featureID] = []
                densityData[assemblyID][chromID][featureID].append((float(position),int(teCount),int(totalTELenPerWindow),float(percentWindow)))
    return(densityData,featureIDs)


def readMethylationData(methylationData):
    methylationDensity = {}
    with open(methylationData,'r') as F:
        for line in F:
            if '#' not in line:
                assemblyID,chromID,position,averageProbability = line.strip().split('\t')
                if assemblyID not in methylationDensity:
                    methylationDensity[assemblyID] = {}
                if chromID not in methylationDensity[assemblyID]:
                    methylationDensity[assemblyID][chromID] = []
                methylationDensity[assemblyID][chromID].append((float(position),float(averageProbability)))
    return(methylationDensity)


def readGCData(gcData):
    gcDataDict = {}
    with open(gcData,'r') as F:
        for line in F:
            if '#' not in line:
                # AH3Mb   chrY    110000000.0     0.5825785116980806      33.59582706777937
                assemblyID,chromID,position,observed_vs_expected,gc_content = line.strip().split('\t')
                if assemblyID not in gcDataDict:
                    gcDataDict[assemblyID] = {}
                if chromID not in gcDataDict[assemblyID]:
                    gcDataDict[assemblyID][chromID] = []
                gcDataDict[assemblyID][chromID].append((float(position),float(observed_vs_expected),float(gc_content)))
    for assemblyID in gcDataDict:
        for chromID in gcDataDict[assemblyID]:
            gcDataDict[assemblyID][chromID].sort(key=lambda x:x[0], reverse=False)
    return(gcDataDict)


def readVariantData(variantData):
    variantDensityDict = {}
    with open(variantData,'r') as F:
        for line in F:
            if 'assemblyID' not in line:
                assemblyID,chromID,position,variantCount,totalVariantLength,totalVariantPercent = line.strip().split('\t')
                if assemblyID not in variantDensityDict:
                    variantDensityDict[assemblyID] = {}
                if chromID not in variantDensityDict[assemblyID]:
                    variantDensityDict[assemblyID][chromID] = []
                variantDensityDict[assemblyID][chromID].append((float(position),int(variantCount),int(totalVariantLength),float(totalVariantPercent)))
    return(variantDensityDict)


def plotData(densityData,methylationDensity,chromosomeLengthDict,colorMap,gcDataDict,geneDensity,lastWindowRemainderDict,featureType,windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict):
    for assemblyID in densityData:
        for chromID in densityData[assemblyID]:
            seqLen = chromosomeLengthDict[assemblyID][chromID]
            lastWindowLength = lastWindowRemainderDict[assemblyID][chromID]
            if assemblyID in dup_variantDensityDict and chromID in dup_variantDensityDict[assemblyID]:
                fig, (ax1, ax6, ax8) = plt.subplots(nrows=3, sharex=True, figsize=(20,10), gridspec_kw=dict(height_ratios=[3, 1, 1]))
            else:
                fig, ax1 = plt.subplots(figsize=(20,10))
            plt.rcParams['font.size'] = 16
            legendDict1 = {}
            legendDict2 = {}
            legendDict3 = {}
            for featureID in densityData[assemblyID][chromID]:
                featureArray = np.array(densityData[assemblyID][chromID][featureID])
                for i in range(len(featureArray)-1):
                    position,teCount,teTotalLen,totalTEPercent = featureArray[i]
                    ax1.bar(position,totalTEPercent,windowSize, align='edge', label=featureID, alpha=0.75, color=colorMap[featureID])
                lastPosition,lastTECount,lastTotalTELen,lastTotalTEPercent = featureArray[-1]
                ax1.bar(lastPosition,lastTotalTEPercent,lastWindowLength, align='edge', label=featureID, alpha=0.75, color=colorMap[featureID])
                ax1.spines['left'].set_color('black')
                ax1.spines['right'].set_color('black')
                ax1.tick_params(axis='y', colors='black')
                ax1.set_ylim(0,100)
                plt.xlabel('Chromosome position (Mb)')
                ax1.set_ylabel('Transposon percent per window (' + str(windowSize) + ' bp)')
                ax1.margins(0)
                legendInfo = mlines.Line2D([], [], color=colorMap[featureID], marker='s',linestyle="None",markersize=10, label=featureID, alpha=1.0, fillstyle='full', markeredgecolor=colorMap[featureID], markeredgewidth=0.0)
                ticks = ax1.get_xticks()/1000000.0
                ax1.set_xticklabels(ticks.astype(int))
                if featureID not in legendDict1:
                    legendDict1[featureID] = legendInfo
            ax2 = ax1.twinx()
            if assemblyID in methylationDensity and chromID in methylationDensity[assemblyID]:
                methylationArray = np.array(methylationDensity[assemblyID][chromID])
                for j in range(len(methylationArray)-2):
                    position,methylationProb = methylationArray[j]
                    # this position modification is to plot the marker in the center of the window
                    position = position + (windowSize/2)
                    nextPosition,nextMethylationProb = methylationArray[j+1]
                    nextPosition = nextPosition + (windowSize/2)
                    ax2.plot((position,nextPosition),(methylationProb,nextMethylationProb), color='black', marker='o', linewidth=4, markersize=4.5)
                # define last position explicitly
                secondToLastPosition,secondToLastMethylationProb = methylationArray[-2]
                lastPosition,lastMethylationProb = methylationArray[-1]
                # this position modification is to plot the marker in the center of the window
                secondToLastPosition = secondToLastPosition + (windowSize/2)
                lastPosition = lastPosition + (lastWindowLength/2)
                ax2.plot((secondToLastPosition,lastPosition),(secondToLastMethylationProb,lastMethylationProb), color='black', marker='o', linewidth=4, markersize=4.5)
                ax2.spines['left'].set_color('black')
                ax2.spines['right'].set_color('black')
                ax2.tick_params(axis='y', colors='black')
                ax2.set_ylabel('Average methylation probability\nper window (' + str(windowSize) + ' bp)')
                ax2.margins(0)
                ax2.set_ylim(0,100)
                # https://stackoverflow.com/questions/10171618/changing-plot-scale-by-a-factor-in-matplotlib
                ticks = ax2.get_xticks()/1000000
                ax2.set_xticklabels(ticks.astype(int))
                legendInfo = mlines.Line2D([], [], color='black', marker='o',markersize=10, label='Methylation probability', alpha=1.0, fillstyle='full', markeredgecolor='black', markeredgewidth=0.0)
                if 'Methylation_probability' not in legendDict1:
                    legendDict1['Methylation_probability'] = legendInfo

            ax3 = ax1.twinx()
            ax3.grid(False)
            ax3.spines.right.set_position(("axes", 1.08))
            if assemblyID in gcDataDict and chromID in gcDataDict[assemblyID]:
                positions,observed_expected,gcContent = zip(*gcDataDict[assemblyID][chromID])
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
                    p3, = ax3.plot((position,nextPosition),(obs_exp,next_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=4.5)
                # explicitly plot last window
                # this position modification is to plot the marker in the center of the window
                secondToLastPos = positionArray[-2]
                lastPos = positionArray[-1]
                secondToLastPos = secondToLastPos + (windowSize/2)
                lastPos = lastPos + (lastWindowLength/2)
                secondToLast_obs_exp = obs_expArray[-2]
                last_obs_exp = obs_expArray[-1]
                ###
                ax3.hlines(y=obs_exp_threshold, xmin=0, xmax=seqLen, linewidth=3, color='#E0A458')
                ###
                ax3.plot((secondToLastPos,lastPos),(secondToLast_obs_exp,last_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=4.5)
                ax3.spines['right'].set_color('#E0A458')
                ax3.tick_params(axis='y', colors='#E0A458')
                ax3.set_ylabel('Observed/expected CpG per window (' + str(windowSize) + ' bp)')
                ax3.yaxis.get_label().set_color(p3.get_color())
                ax3.margins(0)
                ax3.set_ylim(0,1)
                ###
                legendInfo = mlines.Line2D([], [], color='#E0A458', marker='o',markersize=10, label='Observed/expected CpG', alpha=1.0, fillstyle='full', markeredgecolor='#E0A458', markeredgewidth=0.0)
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
                    p4, = ax4.plot((position,nextPosition),(gc,next_gc), color='#824C71', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                secondToLastPos = positionArray[-2]
                lastPos = positionArray[-1]
                secondToLastPos = secondToLastPos + (windowSize/2)
                lastPos = lastPos + (lastWindowLength/2)
                secondToLast_gc = gcContentArray[-2]
                last_gc = gcContentArray[-1]
                ax4.plot((secondToLastPos,lastPos),(secondToLast_gc,last_gc), color='#824C71', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                ax4.spines["left"].set_visible(True)
                ax4.yaxis.set_label_position('left')
                ax4.yaxis.set_ticks_position('left')
                ax4.spines['left'].set_color('#824C71')
                ax4.tick_params(axis='y', colors='#824C71')
                ax4.set_ylabel('GC content (%) per window (' + str(windowSize) + ' bp)')
                ax4.yaxis.get_label().set_color(p4.get_color())
                ax4.margins(0)
                ax4.set_ylim(0,100)
                legendInfo = mlines.Line2D([], [], color='#824C71', marker='o',markersize=10, label='GC content (%)', alpha=1.0, fillstyle='full', markeredgecolor='#824C71', markeredgewidth=0.0, linestyle='dashed')
                if 'GC content (%)' not in legendDict1:
                    legendDict1['GC content (%)'] = legendInfo

            ax5 = ax1.twinx()
            ax5.grid(False)
            ax5.spines.right.set_position(("axes", 1.32))
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
                    p5, = ax5.plot((position,nextPosition),(geneDensityValue,nextGeneDensityValue), color='blue', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                secondToLastPos,secondToLastGeneDensityValue = geneArray[-2]
                secondToLastPos = secondToLastPos + (windowSize/2)
                lastPos,lastGeneDensityValue = geneArray[-1]
                lastPos = lastPos + (lastWindowLength/2)
                ax5.plot((secondToLastPos,lastPos),(secondToLastGeneDensityValue,lastGeneDensityValue), color='blue', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                ax5.spines['right'].set_color('blue')
                ax5.tick_params(axis='y', colors='blue')
                ax5.set_ylabel('Transcript count per window (' + str(windowSize) + ' bp)')
                ax5.yaxis.get_label().set_color(p5.get_color())
                ax5.margins(0)
                legendInfo = mlines.Line2D([], [], color='blue', marker='o',markersize=10, label='Transcript count (' + str(minCount) + '-' + str(maxCount) + ' per 1 Mb)', alpha=1.0, fillstyle='full', markeredgecolor='blue', markeredgewidth=0.0, linestyle='dashed')
                if 'Transcript count' not in legendDict1:
                    legendDict1['Transcript count'] = legendInfo

            variantLabelMap = {'dups':'Duplication', 'transloc':'Translocation', 'inv':'Inversion'}
            variantColorMap = {'dups':'purple', 'transloc':'red', 'inv':'orange'}
            if assemblyID in dup_variantDensityDict and chromID in dup_variantDensityDict[assemblyID]:
                ###
                dupArray = np.array(dup_variantDensityDict[assemblyID][chromID])
                invArray = np.array(inv_variantDensityDict[assemblyID][chromID])
                transArray = np.array(trans_variantDensityDict[assemblyID][chromID])
                for m in range(len(dupArray)-2):
                    dup_position,dup_variantCount,dup_totalVariantLength,dup_totalVariantPercent = dupArray[m]
                    next_dup_position,next_dup_variantCount,next_dup_totalVariantLength,next_dup_totalVariantPercent = dupArray[m+1]
                    dup_position = dup_position + (windowSize/2)
                    next_dup_position = next_dup_position + (windowSize/2)

                    inv_position,inv_variantCount,inv_totalVariantLength,inv_totalVariantPercent = invArray[m]
                    next_inv_position,next_inv_variantCount,next_inv_totalVariantLength,next_inv_totalVariantPercent = invArray[m+1]
                    inv_position = inv_position + (windowSize/2)
                    next_inv_position = next_inv_position + (windowSize/2)

                    trans_position,trans_variantCount,trans_totalVariantLength,trans_totalVariantPercent = transArray[m]
                    next_trans_position,next_trans_variantCount,next_trans_totalVariantLength,next_trans_totalVariantPercent = transArray[m+1]
                    trans_position = trans_position + (windowSize/2)
                    next_trans_position = next_trans_position + (windowSize/2)

                    dup_p6, = ax6.plot((dup_position,next_dup_position),(dup_totalVariantPercent,next_dup_totalVariantPercent), ls='-', color=variantColorMap['dups'], label=variantLabelMap['dups'], linewidth=4)
                    inv_p6, = ax6.plot((inv_position,next_inv_position),(inv_totalVariantPercent,next_inv_totalVariantPercent), ls='-', color=variantColorMap['inv'], label=variantLabelMap['inv'], linewidth=4)
                    trans_p6, = ax6.plot((trans_position,next_trans_position),(trans_totalVariantPercent,next_trans_totalVariantPercent), ls='-', color=variantColorMap['transloc'], label=variantLabelMap['transloc'], linewidth=4)
                ###
                dup_secondToLastPos,dup_secondToLastVariantCount,dup_secondToLastTotalVariantLength,dup_secondToLastTotalVariantPercent = dupArray[-2]
                dup_secondToLastPos = dup_secondToLastPos + (windowSize/2)
                dup_lastPos,dup_lastVariantCount,dup_lastTotalVariantLength,dup_lastTotalVariantPercent = dupArray[-1]
                dup_lastPos = dup_lastPos + (lastWindowLength/2)
                ax6.plot((dup_secondToLastPos,dup_lastPos),(dup_secondToLastTotalVariantPercent,dup_lastTotalVariantPercent), ls='-', color=variantColorMap['dups'], label=variantLabelMap['dups'], linewidth=4)
                ###
                inv_secondToLastPos,inv_secondToLastVariantCount,inv_secondToLastTotalVariantLength,inv_secondToLastTotalVariantPercent = invArray[-2]
                inv_secondToLastPos = inv_secondToLastPos + (windowSize/2)
                inv_lastPos,inv_lastVariantCount,inv_lastTotalVariantLength,inv_lastTotalVariantPercent = invArray[-1]
                inv_lastPos = inv_lastPos + (lastWindowLength/2)
                ax6.plot((inv_secondToLastPos,inv_lastPos),(inv_secondToLastTotalVariantPercent,inv_lastTotalVariantPercent), ls='-', color=variantColorMap['inv'], label=variantLabelMap['inv'], linewidth=4)
                ###
                trans_secondToLastPos,trans_secondToLastVariantCount,trans_secondToLastTotalVariantLength,trans_secondToLastTotalVariantPercent = transArray[-2]
                trans_secondToLastPos = trans_secondToLastPos + (windowSize/2)
                trans_lastPos,trans_lastVariantCount,trans_lastTotalVariantLength,trans_lastTotalVariantPercent = transArray[-1]
                trans_lastPos = trans_lastPos + (lastWindowLength/2)
                ax6.plot((trans_secondToLastPos,trans_lastPos),(trans_secondToLastTotalVariantPercent,trans_lastTotalVariantPercent), ls='-', color=variantColorMap['transloc'], label=variantLabelMap['transloc'], linewidth=4)
                ###
                ax6.set_ylabel('Variant percent\nper window', rotation=45, horizontalalignment='right')
                ax6.set_ylim(0,100)
                dupLegendInfo1 = mlines.Line2D([], [], color=variantColorMap['dups'], marker='o', markersize=10, label=variantLabelMap['dups'] + ' percent', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['dups'], markeredgewidth=0.0)
                if 'Duplication percent' not in legendDict2:
                    legendDict2['Duplication percent'] = dupLegendInfo1
                invLegendInfo1 = mlines.Line2D([], [], color=variantColorMap['inv'], marker='o', markersize=10, label=variantLabelMap['inv'] + ' percent', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['inv'], markeredgewidth=0.0)
                if 'Inversion percent' not in legendDict2:
                    legendDict2['Inversion percent'] = invLegendInfo1
                transLegendInfo1 = mlines.Line2D([], [], color=variantColorMap['transloc'], marker='o', markersize=10, label=variantLabelMap['transloc'] + ' percent', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['transloc'], markeredgewidth=0.0)
                if 'Translocation percent' not in legendDict2:
                    legendDict2['Translocation percent'] = transLegendInfo1
                ax6.margins(0)
            if assemblyID in dup_variantDensityDict and chromID in dup_variantDensityDict[assemblyID]:
                dupArray = np.array(dup_variantDensityDict[assemblyID][chromID])
                invArray = np.array(inv_variantDensityDict[assemblyID][chromID])
                transArray = np.array(trans_variantDensityDict[assemblyID][chromID])
                for n in range(len(dupArray)-2):
                    dup_position,dup_variantCount,dup_totalVariantLength,dup_totalVariantPercent = dupArray[n]
                    next_dup_position,next_dup_variantCount,next_dup_totalVariantLength,next_dup_totalVariantPercent = dupArray[n+1]
                    dup_position = dup_position + (windowSize/2)
                    next_dup_position = next_dup_position + (windowSize/2)
                    
                    inv_position,inv_variantCount,inv_totalVariantLength,inv_totalVariantPercent = invArray[n]
                    next_inv_position,next_inv_variantCount,next_inv_totalVariantLength,next_inv_totalVariantPercent = invArray[n+1]
                    inv_position = inv_position + (windowSize/2)
                    next_inv_position = next_inv_position + (windowSize/2)
                    
                    trans_position,trans_variantCount,trans_totalVariantLength,trans_totalVariantPercent = transArray[n]
                    next_trans_position,next_trans_variantCount,next_trans_totalVariantLength,next_trans_totalVariantPercent = transArray[n+1]
                    trans_position = trans_position + (windowSize/2)
                    next_trans_position = next_trans_position + (windowSize/2)
                    
                    dup_p8 = ax8.plot((dup_position,next_dup_position), (dup_variantCount,next_dup_variantCount), marker='o', linestyle='dashed', color=variantColorMap['dups'], label=variantLabelMap['dups'] + ' count', linewidth=4, alpha=0.5, markersize=4.5)
                    inv_p8 = ax8.plot((inv_position,next_inv_position), (inv_variantCount,next_inv_variantCount), marker='o', linestyle='dashed', color=variantColorMap['inv'], label=variantLabelMap['inv'] + ' count', linewidth=4, alpha=0.5, markersize=4.5)
                    trans_p8 = ax8.plot((trans_position,next_trans_position), (trans_variantCount,next_trans_variantCount), marker='o', linestyle='dashed', color=variantColorMap['transloc'], label=variantLabelMap['transloc'] + ' count', linewidth=4, alpha=0.5, markersize=4.5)
                dup_secondToLastPos,dup_secondToLastVariantCount,dup_secondToLastTotalVariantLength,dup_secondToLastTotalVariantPercent = dupArray[-2]
                dup_secondToLastPos = dup_secondToLastPos + (windowSize/2)
                dup_lastPos,dup_lastVariantCount,dup_lastTotalVariantLength,dup_lastTotalVariantPercent = dupArray[-1]
                dup_lastPos = dup_lastPos + (lastWindowLength/2)
                ax8.plot((dup_secondToLastPos,dup_lastPos),(dup_secondToLastVariantCount,dup_lastVariantCount), ls='-', color=variantColorMap['dups'], label=variantLabelMap['dups'], linewidth=4)
                
                inv_secondToLastPos,inv_secondToLastVariantCount,inv_secondToLastTotalVariantLength,inv_secondToLastTotalVariantPercent = invArray[-2]
                inv_secondToLastPos = inv_secondToLastPos + (windowSize/2)
                inv_lastPos,inv_lastVariantCount,inv_lastTotalVariantLength,inv_lastTotalVariantPercent = invArray[-1]
                inv_lastPos = inv_lastPos + (lastWindowLength/2)
                ax8.plot((inv_secondToLastPos,inv_lastPos),(inv_secondToLastVariantCount,inv_lastVariantCount), ls='-', color=variantColorMap['inv'], label=variantLabelMap['inv'], linewidth=4)

                trans_secondToLastPos,trans_secondToLastVariantCount,trans_secondToLastTotalVariantLength,trans_secondToLastTotalVariantPercent = transArray[-2]
                trans_secondToLastPos = trans_secondToLastPos + (windowSize/2)
                trans_lastPos,trans_lastVariantCount,trans_lastTotalVariantLength,trans_lastTotalVariantPercent = transArray[-1]
                trans_lastPos = trans_lastPos + (lastWindowLength/2)
                ax8.plot((trans_secondToLastPos,trans_lastPos),(trans_secondToLastVariantCount,trans_lastVariantCount), ls='-', color=variantColorMap['transloc'], label=variantLabelMap['transloc'], linewidth=4)

                ax8.set_ylabel('Variant count\nper window', rotation=45, horizontalalignment='right')
                ax8.margins(0)
                dupLegendInfo2 = mlines.Line2D([], [], color=variantColorMap['dups'], marker='o',linestyle="dashed",markersize=10, label=variantLabelMap['dups'] + ' counts', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['dups'], markeredgewidth=0.0)
                if 'Duplication count' not in legendDict3:
                    legendDict3['Duplication count'] = dupLegendInfo2
                invLegendInfo2 = mlines.Line2D([], [], color=variantColorMap['inv'], marker='o',linestyle="dashed",markersize=10, label=variantLabelMap['inv'] + ' counts', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['inv'], markeredgewidth=0.0)
                if 'Inversion count' not in legendDict3:
                    legendDict3['Inversion count'] = invLegendInfo2
                transLegendInfo2 = mlines.Line2D([], [], color=variantColorMap['transloc'], marker='o',linestyle="dashed",markersize=10, label=variantLabelMap['transloc'] + ' counts', alpha=1.0, fillstyle='full', markeredgecolor=variantColorMap['transloc'], markeredgewidth=0.0)
                if 'Translocation count' not in legendDict3:
                    legendDict3['Translocation count'] = transLegendInfo2
            legendList2 = list(legendDict2.values())
            legendList3 = list(legendDict3.values())
            ax6.legend(handles=legendList2, frameon=False, ncol=3, loc='upper center')
            ax8.legend(handles=legendList3, frameon=False, ncol=3, loc='upper center')
            legendList1 = list(legendDict1.values())
            ax1.legend(handles=legendList1, frameon=False, ncol=2, loc='upper center')
            ax1.set_title(assemblyID + " " + chromID + " " + str(seqLen) + " bp")
            plt.savefig(assemblyID + "_" + chromID + "_" + featureType + "_distribution.png",dpi=600)
            plt.savefig(assemblyID + "_" + chromID + "_" + featureType + "_distribution.svg",)
            plt.close()



def plotTimes(compiledTimes,methylationDensity,chromosomeLengthDict,colorMap,gcDataDict,geneDensity,lastWindowRemainderDict,featureType,windowSize,obs_exp_threshold):
    colorDict = {'0_5':'#AA3377', '6_10':'#EE6677', '11_15':'#66CCEE', '16_20':'#228833', '21_25':'#4477AA'}
    labelDict = {'0_5':'<6 mya', '6_10':'>=6 to <11 mya', '11_15':'>=11 to <16 mya', '16_20':'>=16 to <21 mya', '21_25':'>=21 to <26 mya'}
    for assemblyID in compiledTimes:
        for chromID in compiledTimes[assemblyID]:
            seqLen = chromosomeLengthDict[assemblyID][chromID]
            lastWindowLength = lastWindowRemainderDict[assemblyID][chromID]
            for featureID in compiledTimes[assemblyID][chromID]:
                # fig, ax1 = plt.subplots(figsize=(25,12.5))
                fig, ax1 = plt.subplots(figsize=(20,10))
                plt.rcParams['font.size'] = 16
                legendDict = {}
                for timeID in compiledTimes[assemblyID][chromID][featureID]:
                    featureArray = np.array(compiledTimes[assemblyID][chromID][featureID][timeID])
                    for i in range(len(featureArray)-1):
                        position,density,totalTELen,totalTEPercent,avgTime,specificTimePercent = featureArray[i]
                        ax1.bar(position,specificTimePercent,windowSize, align='edge', label=featureID, alpha=0.75, color=colorDict[timeID])
                    lastPosition,lastDensity,lastTotalTELen,lastTotalTEPercent,lastAvgTime,lastSpecificTimePercent = featureArray[-1]
                    ax1.bar(lastPosition,lastSpecificTimePercent,lastWindowLength, align='edge', label=featureID, alpha=0.75, color=colorDict[timeID])
                    ax1.spines['left'].set_color('black')
                    ax1.spines['right'].set_color('black')
                    ax1.tick_params(axis='y', colors='black')
                    plt.xlabel('Chromosome position (Mb)')
                    plt.ylabel('Transposon percent per window (' + str(windowSize) + ' bp)')
                    #plt.margins(0)
                    ax1.margins(0)
                    ax1.set_ylim(0,100)
                    legendInfo = mlines.Line2D([], [], color=colorDict[timeID], marker='s',linestyle="None",markersize=10, label=labelDict[timeID], alpha=1.0, fillstyle='full', markeredgecolor=colorDict[timeID], markeredgewidth=0.0)
                    ticks = ax1.get_xticks()/1000000.0
                    ax1.set_xticklabels(ticks.astype(int))
                    if timeID not in legendDict:
                        legendDict[timeID] = legendInfo

                ax2 = ax1.twinx()
                if assemblyID in methylationDensity and chromID in methylationDensity[assemblyID]:
                    methylationArray = np.array(methylationDensity[assemblyID][chromID])
                    for j in range(len(methylationArray)-2):
                        position,methylationProb = methylationArray[j]
                        # this position modification is to plot the marker in the center of the window
                        position = position + (windowSize/2)
                        nextPosition,nextMethylationProb = methylationArray[j+1]
                        nextPosition = nextPosition + (windowSize/2)
                        ax2.plot((position,nextPosition),(methylationProb,nextMethylationProb), color='black', marker='o', linewidth=4, markersize=4.5)
                    # define last position explicitly
                    secondToLastPosition,secondToLastMethylationProb = methylationArray[-2]
                    lastPosition,lastMethylationProb = methylationArray[-1]
                    # this position modification is to plot the marker in the center of the window
                    secondToLastPosition = secondToLastPosition + (windowSize/2)
                    lastPosition = lastPosition + (lastWindowLength/2)
                    ax2.plot((secondToLastPosition,lastPosition),(secondToLastMethylationProb,lastMethylationProb), color='black', marker='o', linewidth=4, markersize=4.5)
                    ax2.spines['left'].set_color('black')
                    ax2.spines['right'].set_color('black')
                    ax2.tick_params(axis='y', colors='black')
                    ax2.set_ylabel('Average methylation probability\nper window (' + str(windowSize) + ' bp)')
                    #plt.margins(0)
                    ax2.margins(0)
                    ax2.set_ylim(0,100)
                    # https://stackoverflow.com/questions/10171618/changing-plot-scale-by-a-factor-in-matplotlib
                    ticks = ax2.get_xticks()/1000000
                    ax2.set_xticklabels(ticks.astype(int))
                    legendInfo = mlines.Line2D([], [], color='black', marker='o',markersize=10, label='Methylation probability', alpha=1.0, fillstyle='full', markeredgecolor='black', markeredgewidth=0.0)
                    if 'Methylation_probability' not in legendDict:
                        legendDict['Methylation_probability'] = legendInfo

                ax3 = ax1.twinx()
                ax3.grid(False)
                ax3.spines.right.set_position(("axes", 1.08))
                if assemblyID in gcDataDict and chromID in gcDataDict[assemblyID]:
                    positions,observed_expected,gcContent = zip(*gcDataDict[assemblyID][chromID])
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
                        p3, = ax3.plot((position,nextPosition),(obs_exp,next_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=4.5)
                    # explicitly plot last window
                    # this position modification is to plot the marker in the center of the window
                    secondToLastPos = positionArray[-2]
                    lastPos = positionArray[-1]
                    secondToLastPos = secondToLastPos + (windowSize/2)
                    lastPos = lastPos + (lastWindowLength/2)
                    secondToLast_obs_exp = obs_expArray[-2]
                    last_obs_exp = obs_expArray[-1]
                    ax3.hlines(y=obs_exp_threshold, xmin=0, xmax=seqLen, linewidth=3, color='#E0A458')
                    ax3.plot((secondToLastPos,lastPos),(secondToLast_obs_exp,last_obs_exp), color='#E0A458', marker='o', linewidth=4, markersize=4.5)
                    ax3.spines['right'].set_color('#E0A458')
                    ax3.tick_params(axis='y', colors='#E0A458')
                    ax3.set_ylabel('Observed/expected CpG per window (' + str(windowSize) + ' bp)')
                    ax3.yaxis.get_label().set_color(p3.get_color())
                    #plt.margins(0)
                    ax3.margins(0)
                    ax3.set_ylim(0,1)
                    ###
                    legendInfo = mlines.Line2D([], [], color='#E0A458', marker='o',markersize=10, label='Observed/expected CpG', alpha=1.0, fillstyle='full', markeredgecolor='#E0A458', markeredgewidth=0.0)
                    if 'Observed/expected CpG' not in legendDict:
                        legendDict['Observed/expected CpG'] = legendInfo
                        
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
                        p4, = ax4.plot((position,nextPosition),(gc,next_gc), color='#824C71', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                    secondToLastPos = positionArray[-2]
                    lastPos = positionArray[-1]
                    secondToLastPos = secondToLastPos + (windowSize/2)
                    lastPos = lastPos + (lastWindowLength/2)
                    secondToLast_gc = gcContentArray[-2]
                    last_gc = gcContentArray[-1]
                    ax4.plot((secondToLastPos,lastPos),(secondToLast_gc,last_gc), color='#824C71', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                    ax4.spines["left"].set_visible(True)
                    ax4.yaxis.set_label_position('left')
                    ax4.yaxis.set_ticks_position('left')
                    ax4.spines['left'].set_color('#824C71')
                    ax4.tick_params(axis='y', colors='#824C71')
                    ax4.set_ylabel('GC content (%) per window (' + str(windowSize) + ' bp)')
                    ax4.yaxis.get_label().set_color(p4.get_color())
                    #plt.margins(0)
                    ax4.margins(0)
                    ax4.set_ylim(0,100)
                    legendInfo = mlines.Line2D([], [], color='#824C71', marker='o',markersize=10, label='GC content (%)', alpha=1.0, fillstyle='full', markeredgecolor='#824C71', markeredgewidth=0.0, linestyle='dashed')
                    if 'GC content (%)' not in legendDict:
                        legendDict['GC content (%)'] = legendInfo
                
                ax5 = ax1.twinx()
                ax5.grid(False)
                ax5.spines.right.set_position(("axes", 1.32))
                if assemblyID in geneDensity and chromID in geneDensity[assemblyID]:
                    positions,counts = zip(*geneDensity[assemblyID][chromID])
                    minCount = min(counts)
                    maxCount = max(counts)
                    geneArray = np.array(geneDensity[assemblyID][chromID])
                    for k in range(len(geneArray)-2):
                        position,geneCount = geneArray[k]
                        nextPosition,nextGeneCount = geneArray[k+1]
                        position = position + (windowSize/2)
                        nextPosition = nextPosition + (windowSize/2)
                        p5, = ax5.plot((position,nextPosition),(geneCount,nextGeneCount), color='blue', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                    secondToLastPos,secondToLastGeneCount = geneArray[-2]
                    secondToLastPos = secondToLastPos + (windowSize/2)
                    lastPos,lastGeneCount = geneArray[-1]
                    lastPos = lastPos + (lastWindowLength/2)
                    ax5.plot((secondToLastPos,lastPos),(secondToLastGeneCount,lastGeneCount), color='blue', marker='o', linewidth=4, markersize=4.5, linestyle='dashed')
                    ax5.spines['right'].set_color('blue')
                    ax5.tick_params(axis='y', colors='blue')
                    ax5.set_ylabel('Transcript count per window (' + str(windowSize) + ' bp)')
                    ax5.yaxis.get_label().set_color(p5.get_color())
                    #plt.margins(0)
                    ax5.margins(0)
                    legendInfo = mlines.Line2D([], [], color='blue', marker='o',markersize=10, label='Transcript count (' + str(minCount) + '-' + str(maxCount) + ' per 1 Mb)', alpha=1.0, fillstyle='full', markeredgecolor='blue', markeredgewidth=0.0, linestyle='dashed')
                    if 'Transcript count' not in legendDict:
                        legendDict['Transcript count'] = legendInfo

                legendList = list(legendDict.values())
                plt.title(assemblyID + " " + chromID + " " + featureID + " " + str(seqLen) + " bp")
                plt.legend(handles=legendList, frameon=False, ncol=2, loc='upper center')
                plt.savefig(assemblyID + "_" + chromID + "_" + featureID + "_subset_time_densities.png",dpi=600)
                plt.savefig(assemblyID + "_" + chromID + "_" + featureID + "_subset_time_densities.svg")
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
                
            
usage = "Usage: " + sys.argv[0] + " <gene data> <DNA transposon density/count data> <DNA transposon time data> <long terminal retrotransposon density/count data> <long terminal retrotransposon time data> <methylation data> <GC data> <duplication variant data> <inversion variant data> <translocation variant data> <chromosome lengths> <assembly ID> <window size, e.g. 1000000>\n"
if len(sys.argv) != 14:
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
assemblyID = sys.argv[12]
windowSize = sys.argv[13]

windowSize = int(windowSize)
obs_exp_threshold = 0.6

# chromosome lengths
chromosomeLengthDict = readChromosomeLengthData(chromosomeLengths)
lastWindowRemainderDict = calculateLastWindowLength(chromosomeLengthDict,windowSize)

# gene counts
geneDensity = readGeneData(geneData)

# transposable elements densities/counts
teDensityData,teFeatureIDs = readTEDensityData(teData)
ltrDensityData,ltrFeatureIDs = readTEDensityData(ltrData)

# transposable element times
teTimeDataDict = readTETimeData(teTimeData)
ltrTimeDataDict = readTETimeData(ltrTimeData)

# gc data
gcDataDict = readGCData(gcData)

# methylation data
methylationDensity = readMethylationData(methylationData)

# structural variants
dup_variantDensityDict = readVariantData(dupData)
inv_variantDensityDict = readVariantData(invData)
trans_variantDensityDict = readVariantData(translocData)

# color maps for plotting
teColorMap = createColorMapDict(teFeatureIDs,'TEs')
ltrColorMap = createColorMapDict(ltrFeatureIDs,'LTRs')

# plot data
plotTimes(teTimeDataDict,methylationDensity,chromosomeLengthDict,teColorMap,gcDataDict,geneDensity,lastWindowRemainderDict,'TEs',windowSize,obs_exp_threshold)
plotData(teDensityData,methylationDensity,chromosomeLengthDict,teColorMap,gcDataDict,geneDensity,lastWindowRemainderDict,'TEs',windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict)

plotTimes(ltrTimeDataDict,methylationDensity,chromosomeLengthDict,ltrColorMap,gcDataDict,geneDensity,lastWindowRemainderDict,'LTRs',windowSize,obs_exp_threshold)
plotData(ltrDensityData,methylationDensity,chromosomeLengthDict,ltrColorMap,gcDataDict,geneDensity,lastWindowRemainderDict,'LTRs',windowSize,obs_exp_threshold,dup_variantDensityDict,inv_variantDensityDict,trans_variantDensityDict)
