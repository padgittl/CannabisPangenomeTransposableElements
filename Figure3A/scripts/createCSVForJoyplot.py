import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import sem
from joypy import joyplot
import pandas as pd
import random
import seaborn as sns
from matplotlib import cm


########
# MAIN #
########


def readChemotypeFile(chemotypeFile):
    populationDict = {}
    # high-cannabinoid hemp
    cleanPopIDs = {'mj':'mj', 'hc_hemp':'hc hemp',
                   'feral':'feral', 'hemp':'hemp', 'asian_hemp':'Asian hemp', 'F1':'F1'}
    with open(chemotypeFile,'r') as F:
        for line in F:
            if 'Chemotype' not in line and 'EXCLUDED' not in line and ',,' not in line:
                assemblyID,chemotypeID,populationID = line.strip().split(',')
                if populationID in cleanPopIDs:
                    newPopID = cleanPopIDs[populationID]
                else:
                    print(line)
                    print(populationID,"not in cleanPopIDs")
                    sys.exit()
                if assemblyID not in populationDict:
                    populationDict[assemblyID] = (chemotypeID,newPopID)
    return(populationDict)


def readFileList(fileList):
    dataList = []
    with open(fileList,'r') as F:
        for line in F:
            fileName = line.strip()
            dataList.append(fileName)
    return(dataList)

'''
Count        bpMasked    %masked
TIR                    --           --           --   
    CACTA              97530        27175721     3.51% 
    Mutator            63593        37094405     4.79% 
    PIF_Harbinger      15127        10493347     1.35% 
    Tc1_Mariner        4810         2015085      0.26% 
    hAT                56173        14872409     1.92% 
nonTIR                 --           --           --   
    helitron           87485        21694748     2.80% 
teColors = {'CACTA_TIR_transposon':'#000000','Mutator_TIR_transposon':'#E69F00','PIF_Harbinger_TIR_transposon':'#0072B2','Tc1_Mariner_TIR_transposon':'#009E73','hAT_TIR_transposon':'#F0E442','helitron':'#CC79A7'}
    ltrColors = {'Ty1_LTR_retrotransposon':'#542788', 'unknown_LTR_retrotransposon':'#998ec3', 'Ty3_LTR_retrotransposon':'#f1a340'}
Copia              153448       125160457    16.16% 
Gypsy              130032       157346446    20.32% 
unknown            272589       127529528    16.47% 
'''
def readSumFiles(sumFileDataList,populationDict):
    CACTAList = []
    MutatorList = []
    PIF_HarbingerList = []
    Tc1_MarinerList = []
    hATList = []
    helitronList = []
    Ty1List = []
    Ty3List = []
    unknownLTRList = []
    totalPercentList = []
    dataDict = {}
    total_dictForPrinting = {}
    dna_dictForPrinting = {}
    ltr_dictForPrinting = {}
    for fileID in sumFileDataList:
        splitFileID = os.path.basename(fileID) 
        assemblyID,suffix = splitFileID.split('.unmasked.fasta.mod.EDTA.TEanno.sum')
        if assemblyID in populationDict:
            chemotypeID,populationID = populationDict[assemblyID]
            if populationID not in dataDict:
                dataDict[populationID] = {}
                
            if populationID not in total_dictForPrinting:
                total_dictForPrinting[populationID] = {}
            if assemblyID not in total_dictForPrinting[populationID]:
                total_dictForPrinting[populationID][assemblyID] = []

            if populationID not in dna_dictForPrinting:
                dna_dictForPrinting[populationID] = {}
            if assemblyID not in dna_dictForPrinting[populationID]:
                dna_dictForPrinting[populationID][assemblyID] = []

            if populationID not in ltr_dictForPrinting:
                ltr_dictForPrinting[populationID] = {}
            if assemblyID not in ltr_dictForPrinting[populationID]:
                ltr_dictForPrinting[populationID][assemblyID] = []
            with open(fileID,'r') as F:
                for line in F:
                    # Total                  938538       544067936    69.72%
                    if 'Total' in line and '%' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        totalPercentList.append(percentMasked)
                        if 'Total' not in dataDict[populationID]:
                            dataDict[populationID]['Total'] = []
                        dataDict[populationID]['Total'].append((assemblyID,percentMasked))
                        total_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    if 'Copia' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        #print(featureID)
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        Ty1List.append(percentMasked)
                        if 'Ty1' not in dataDict[populationID]:
                            dataDict[populationID]['Ty1'] = []
                        dataDict[populationID]['Ty1'].append((assemblyID,percentMasked))
                        ltr_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                        ###
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    if 'Gypsy' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        Ty3List.append(percentMasked)
                        if 'Ty3' not in dataDict[populationID]:
                            dataDict[populationID]['Ty3'] = []
                        dataDict[populationID]['Ty3'].append((assemblyID,percentMasked))
                        ltr_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                        ###
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    if 'unknown' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        unknownLTRList.append(percentMasked)
                        if 'unknown' not in dataDict[populationID]:
                            dataDict[populationID]['unknown'] = []
                        dataDict[populationID]['unknown'].append((assemblyID,percentMasked))
                        ltr_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                        ###
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    if 'CACTA' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        CACTAList.append(percentMasked)
                        if 'CACTA' not in dataDict[populationID]:
                            dataDict[populationID]['CACTA'] = []
                        dataDict[populationID]['CACTA'].append((assemblyID,percentMasked))
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    elif 'Mutator' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        MutatorList.append(percentMasked)
                        if 'Mutator' not in dataDict[populationID]:
                            dataDict[populationID]['Mutator'] = []
                        dataDict[populationID]['Mutator'].append((assemblyID,percentMasked))
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    elif 'PIF_Harbinger' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        PIF_HarbingerList.append(percentMasked)
                        if 'PIF_Harbinger' not in dataDict[populationID]:
                            dataDict[populationID]['PIF_Harbinger'] = []
                        dataDict[populationID]['PIF_Harbinger'].append((assemblyID,percentMasked))
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    elif 'Tc1_Mariner' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        Tc1_MarinerList.append(percentMasked)
                        if 'Tc1_Mariner' not in dataDict[populationID]:
                            dataDict[populationID]['Tc1_Mariner'] = []
                        dataDict[populationID]['Tc1_Mariner'].append((assemblyID,percentMasked))
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    elif 'hAT' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        hATList.append(percentMasked)
                        if 'hAT' not in dataDict[populationID]:
                            dataDict[populationID]['hAT'] = []
                        dataDict[populationID]['hAT'].append((assemblyID,percentMasked))
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
                    elif 'helitron' in line:
                        featureID,featureCount,bpMasked,percentMasked = line.strip().split()
                        percentMasked,extra = percentMasked.split('%')
                        percentMasked = float(percentMasked)
                        helitronList.append(percentMasked)
                        if 'helitron' not in dataDict[populationID]:
                            dataDict[populationID]['helitron'] = []
                        dataDict[populationID]['helitron'].append((assemblyID,percentMasked))
                        dna_dictForPrinting[populationID][assemblyID].append((featureID,percentMasked))
    avgTotalRepeatPercent = sum(totalPercentList) / len(totalPercentList)
    avgTotalRepeatPercent = round(avgTotalRepeatPercent,2)
    totalStdDev = np.std(totalPercentList)
    totalStdDev = round(totalStdDev,2)
    totalSEM = sem(totalPercentList)
    totalSEM = round(totalSEM,2)
    #print("Average total repeat %: ",avgTotalRepeatPercent," SD=",totalStdDev," SEM=",totalSEM)
    
    avgCACTARepeatPercent = sum(CACTAList) / len(CACTAList)
    avgCACTARepeatPercent = round(avgCACTARepeatPercent,2)
    CACTAStdDev = np.std(CACTAList)
    CACTAStdDev = round(CACTAStdDev,2)
    cactaSEM = sem(CACTAList)
    cactaSEM = round(cactaSEM,2)
    #print("Average CACTA repeat %: ",avgCACTARepeatPercent," SD=",CACTAStdDev," SEM=",cactaSEM)

    avgMutatorRepeatPercent = sum(MutatorList) / len(MutatorList)
    avgMutatorRepeatPercent = round(avgMutatorRepeatPercent,2)
    MutatorStdDev = np.std(MutatorList)
    MutatorStdDev = round(MutatorStdDev,2)
    mutatorSEM = sem(MutatorList)
    mutatorSEM = round(mutatorSEM,2)
    #print("Average Mutator repeat %: ",avgMutatorRepeatPercent," SD=",MutatorStdDev," SEM=",mutatorSEM)

    avgPIF_HarbingerRepeatPercent = sum(PIF_HarbingerList) / len(PIF_HarbingerList)
    avgPIF_HarbingerRepeatPercent = round(avgPIF_HarbingerRepeatPercent,2)
    PIF_HarbingerStdDev = np.std(PIF_HarbingerList)
    PIF_HarbingerStdDev = round(PIF_HarbingerStdDev,2)
    harbingerSEM = sem(PIF_HarbingerList)
    harbingerSEM = round(harbingerSEM,2)
    #print("Average PIF_Harbinger repeat %: ",avgPIF_HarbingerRepeatPercent," SD=",PIF_HarbingerStdDev," SEM=",harbingerSEM)

    avgTc1_MarinerRepeatPercent = sum(Tc1_MarinerList) / len(Tc1_MarinerList)
    avgTc1_MarinerRepeatPercent = round(avgTc1_MarinerRepeatPercent,2)
    Tc1_MarinerStdDev = np.std(Tc1_MarinerList)
    Tc1_MarinerStdDev = round(Tc1_MarinerStdDev,2)
    marinerSEM = sem(Tc1_MarinerList)
    marinerSEM = round(marinerSEM,2)
    #print("Average Tc1_Mariner repeat %: ",avgTc1_MarinerRepeatPercent," SD=",Tc1_MarinerStdDev," SEM=",marinerSEM)

    avghATRepeatPercent = sum(hATList) / len(hATList)
    avghATRepeatPercent = round(avghATRepeatPercent,2)
    hATStdDev = np.std(hATList)
    hATStdDev = round(hATStdDev,2)
    hAT_SEM = sem(hATList)
    hAT_SEM = round(hAT_SEM,2)
    #print("Average hAT repeat %: ",avghATRepeatPercent," SD=",hATStdDev," SEM=",hAT_SEM)

    avghelitronRepeatPercent = sum(helitronList) / len(helitronList)
    avghelitronRepeatPercent = round(avghelitronRepeatPercent,2)
    helitronStdDev = np.std(helitronList)
    helitronStdDev = round(helitronStdDev,2)
    helitronSEM = sem(helitronList)
    helitronSEM = round(helitronSEM,2)
    #print("Average helitron repeat %: ",avghelitronRepeatPercent," SD=",helitronStdDev," SEM=",helitronSEM)

    avgTy1RepeatPercent = sum(Ty1List) / len(Ty1List)
    avgTy1RepeatPercent = round(avgTy1RepeatPercent,2)
    Ty1StdDev = np.std(Ty1List)
    Ty1StdDev = round(Ty1StdDev,2)
    ty1SEM = sem(Ty1List)
    ty1SEM = round(ty1SEM,2)
    #print("Average Ty1 repeat %: ",avgTy1RepeatPercent," SD=",Ty1StdDev," SEM=",ty1SEM)

    avgTy3RepeatPercent = sum(Ty3List) / len(Ty3List)
    avgTy3RepeatPercent = round(avgTy3RepeatPercent,2)
    Ty3StdDev = np.std(Ty3List)
    Ty3StdDev = round(Ty3StdDev,2)
    ty3SEM = sem(Ty3List)
    ty3SEM = round(ty3SEM,2)
    #print("Average Ty3 repeat %: ",avgTy3RepeatPercent," SD=",Ty3StdDev," SEM=",ty3SEM)

    avgunknownLTRRepeatPercent = sum(unknownLTRList) / len(unknownLTRList)
    avgunknownLTRRepeatPercent = round(avgunknownLTRRepeatPercent,2)
    unknownLTRStdDev = np.std(unknownLTRList)
    unknownLTRStdDev = round(unknownLTRStdDev,2)
    unknownLTR_SEM = sem(unknownLTRList)
    unknownLTR_SEM = round(unknownLTR_SEM,2)
    #print("Average unknownLTR repeat %: ",avgunknownLTRRepeatPercent," SD=",unknownLTRStdDev," SEM=",unknownLTR_SEM)
    return(dna_dictForPrinting,ltr_dictForPrinting,total_dictForPrinting)


def prepareOutputFile(dictForPrinting,keyword):
    setList = []
    finalData = {}
    OUT = open(keyword + '_percentMasked.csv','w')
    '''
    mj
    hc_hemp
    feral
    hemp
    asian_hemp
    F1
    '''
    for populationID in dictForPrinting:
        for assemblyID in dictForPrinting[populationID]:
            featureIDs,percentMasked = zip(*dictForPrinting[populationID][assemblyID])
            joinedPercentMasked = ",".join(map(str, percentMasked))
            joinedPercentMasked = joinedPercentMasked + ',' + populationID
            setList.append(featureIDs)
            # finalData.append(joinedPercentMasked)
            # .sort(key = lambda x:x[0], reverse=False)
            if populationID not in finalData:
                finalData[populationID] = []
            finalData[populationID].append(joinedPercentMasked)
            
    dataSet = set(setList)
    if len(dataSet) == 1:
        dataSetList = list(dataSet)
        for i in dataSetList:
            columns = ','.join(i)
            columns = columns.replace('Copia','Ty1-LTR')
            columns = columns.replace('Gypsy','Ty3-LTR')
            columns = columns.replace('unknown','unknown-LTR')
            columns = columns + ',Names'
            # finalData.append(columns)
            OUT.write("%s\n" % (columns))
    else:
        print("TE IDs are out of order")
        sys.exit()
    #finalDataList = finalData[::-1]
    for popID in finalData:
        for i in finalData[popID]:
            OUT.write("%s\n" % (i))
    



usage = "Usage: " + sys.argv[0] + " <sum file list> <chemotype-population ID file>\n"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

sumFileList = sys.argv[1]
chemotypeFile = sys.argv[2]

populationDict = readChemotypeFile(chemotypeFile)

sumFileDataList = readFileList(sumFileList)
dna_dictForPrinting,ltr_dictForPrinting,total_dictForPrinting = readSumFiles(sumFileDataList,populationDict)

prepareOutputFile(dna_dictForPrinting,'dnaTE')
prepareOutputFile(ltr_dictForPrinting,'ltr')
prepareOutputFile(total_dictForPrinting,'total')
