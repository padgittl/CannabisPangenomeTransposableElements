import sys, re, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

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


# https://github.com/oushujun/EDTA/issues/233
# "You may grep LTR_retrotransposon *.fasta.mod.EDTA.intact.gff3|less -S to get only intact LTRs."
# "For intact LTR information in EDTA results, you can find the identity info from the last column of the gff3 file. You can use this to calculate the age yourself, with T = K / 2µ = (1 - identity) / 2µ, where µ is the mutation rate of your species in the unit of per bp per year. By default, LTR_retriever uses µ = 1.3e-8 from rice."
def readTEAnnoGFF(teAnnoGFFDataList):
    lambda_value = 6.1e-09
    # distanceThreshold = 0.9987 # 100k years ago
    # >>> (float(1-0.9987799)/(2*6.1e-09))
    # 100008.19672130706
    # distanceThreshold = 0.987
    ## distanceThreshold = 0.9987799 # ~100 kya
    # distanceThreshold = 0.9986 # ~114 kya
    distanceThreshold = 0.998 # ~150 kya
    # distanceThreshold = 0.98 # ~1.6 mya
    # distanceThreshold = 0.9935 # ~500 kya
    # distanceThreshold = 0.997 # ~250 kya
    thousandsTimeDivisor = 1000.0
    millionsTimeDivisor  = 1000000.0
    maxTime = ((1-distanceThreshold) / (2*lambda_value)) / thousandsTimeDivisor
    # print(maxTime)
    
    ty1TimeDictFullNonIntact = {}
    ty3TimeDictFullNonIntact = {}
    unknownTimeDictFullNonIntact = {}
    CACTATimeDictFullNonIntact = {}
    MutatorTimeDictFullNonIntact = {}
    HarbingerTimeDictFullNonIntact = {}
    MarinerTimeDictFullNonIntact = {}
    hATTimeDictFullNonIntact = {}
    helitronTimeDictFullNonIntact = {}

    ty1TimeDictSubsetNonIntact = {}
    ty3TimeDictSubsetNonIntact = {}
    unknownTimeDictSubsetNonIntact = {}
    CACTATimeDictSubsetNonIntact = {}
    MutatorTimeDictSubsetNonIntact = {}
    HarbingerTimeDictSubsetNonIntact = {}
    MarinerTimeDictSubsetNonIntact = {}
    hATTimeDictSubsetNonIntact = {}
    helitronTimeDictSubsetNonIntact = {}

    ty1TimeDictFullIntact = {}
    ty3TimeDictFullIntact = {}
    unknownTimeDictFullIntact = {}
    CACTATimeDictFullIntact = {}
    MutatorTimeDictFullIntact = {}
    HarbingerTimeDictFullIntact = {}
    MarinerTimeDictFullIntact = {}
    hATTimeDictFullIntact = {}
    helitronTimeDictFullIntact = {}

    ty1TimeDictSubsetIntact = {}
    ty3TimeDictSubsetIntact = {}
    unknownTimeDictSubsetIntact = {}
    CACTATimeDictSubsetIntact = {}
    MutatorTimeDictSubsetIntact = {}
    HarbingerTimeDictSubsetIntact = {}
    MarinerTimeDictSubsetIntact = {}
    hATTimeDictSubsetIntact = {}
    helitronTimeDictSubsetIntact = {}

    OUT_full = open('fullSet_TEs_times.txt','w')
    OUT_subset = open('subset_TEs_times.txt','w')
    for gffFile in teAnnoGFFDataList:
        with open(gffFile,'r') as F:
            for line in F:
                if not line.startswith('#'):
                    scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                    if 'chr' in scaffoldID:
                        assemblyID,chrID = scaffoldID.split('.')
                        if feature == 'Copia_LTR_retrotransposon':
                            # Identity=0.903;Method=homology
                            # E_00003284;Classification=LTR/Copia;Sequence_ontology=SO:0002264;ltr_identity=0.9802;Method=structural;
                            # suffix,ltrIdentity = newAttribute[5].split('ltr_identity=')
                            if 'Method=homology' in attribute:
                                newAttribute = attribute.split(';')
                                suffix,identity = newAttribute[4].split('Identity=')
                                identity = identity.strip()
                                if identity != 'NA':
                                    identity = float(identity)
                                    dist = 1-identity
                                    time = (float(dist) / (2*lambda_value))
                                    time_in_thousands = float(time) / thousandsTimeDivisor
                                    time_in_millions = float(time) / millionsTimeDivisor
                                    if 'Copia' not in ty1TimeDictFullNonIntact:
                                        ty1TimeDictFullNonIntact['Copia'] = {}
                                    if chrID not in ty1TimeDictFullNonIntact['Copia']:
                                        ty1TimeDictFullNonIntact['Copia'][chrID] = {}
                                    if assemblyID not in ty1TimeDictFullNonIntact['Copia'][chrID]:
                                        ty1TimeDictFullNonIntact['Copia'][chrID][assemblyID] = []
                                    ty1TimeDictFullNonIntact['Copia'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Copia' not in ty1TimeDictSubsetNonIntact:
                                            ty1TimeDictSubsetNonIntact['Copia'] = {}
                                        if chrID not in ty1TimeDictSubsetNonIntact['Copia']:
                                            ty1TimeDictSubsetNonIntact['Copia'][chrID] = {}
                                        if assemblyID not in ty1TimeDictSubsetNonIntact['Copia'][chrID]:
                                            ty1TimeDictSubsetNonIntact['Copia'][chrID][assemblyID] = []
                                        ty1TimeDictSubsetNonIntact['Copia'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                            if 'Method=structural' in attribute:
                                newAttribute = attribute.split(';')
                                suffix,identity = newAttribute[5].split('ltr_identity=')
                                identity = identity.strip()
                                if identity != 'NA':
                                    identity = float(identity)
                                    dist = 1-identity
                                    time = (float(dist) / (2*lambda_value))
                                    time_in_thousands = float(time) / thousandsTimeDivisor
                                    time_in_millions = float(time) / millionsTimeDivisor
                                    if 'Copia' not in ty1TimeDictFullIntact:
                                        ty1TimeDictFullIntact['Copia'] = {}
                                    if chrID not in ty1TimeDictFullIntact['Copia']:
                                        ty1TimeDictFullIntact['Copia'][chrID] = {}
                                    if assemblyID not in ty1TimeDictFullIntact['Copia'][chrID]:
                                        ty1TimeDictFullIntact['Copia'][chrID][assemblyID] = []
                                    ty1TimeDictFullIntact['Copia'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Copia' not in ty1TimeDictSubsetIntact:
                                            ty1TimeDictSubsetIntact['Copia'] = {}
                                        if chrID not in ty1TimeDictSubsetIntact['Copia']:
                                            ty1TimeDictSubsetIntact['Copia'][chrID] = {}
                                        if assemblyID not in ty1TimeDictSubsetIntact['Copia'][chrID]:
                                            ty1TimeDictSubsetIntact['Copia'][chrID][assemblyID] = []
                                        ty1TimeDictSubsetIntact['Copia'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if feature == 'Gypsy_LTR_retrotransposon':
                            if 'Method=homology' in attribute:
                                newAttribute = attribute.split(';')
                                suffix,identity = newAttribute[4].split('Identity=')
                                identity = identity.strip()
                                if identity != 'NA':
                                    identity = float(identity)
                                    dist = 1-identity
                                    time = (float(dist) / (2*lambda_value))
                                    time_in_thousands = float(time) / thousandsTimeDivisor
                                    time_in_millions = float(time) / millionsTimeDivisor
                                    if 'Gypsy' not in ty3TimeDictFullNonIntact:
                                        ty3TimeDictFullNonIntact['Gypsy'] = {}
                                    if chrID not in ty3TimeDictFullNonIntact['Gypsy']:
                                        ty3TimeDictFullNonIntact['Gypsy'][chrID] = {}
                                    if assemblyID not in ty3TimeDictFullNonIntact['Gypsy'][chrID]:
                                        ty3TimeDictFullNonIntact['Gypsy'][chrID][assemblyID] = []
                                    ty3TimeDictFullNonIntact['Gypsy'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Gypsy' not in ty3TimeDictSubsetNonIntact:
                                            ty3TimeDictSubsetNonIntact['Gypsy'] = {}
                                        if chrID not in ty3TimeDictSubsetNonIntact['Gypsy']:
                                            ty3TimeDictSubsetNonIntact['Gypsy'][chrID] = {}
                                        if assemblyID not in ty3TimeDictSubsetNonIntact['Gypsy'][chrID]:
                                            ty3TimeDictSubsetNonIntact['Gypsy'][chrID][assemblyID] = []
                                        ty3TimeDictSubsetNonIntact['Gypsy'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                            if 'Method=structural' in attribute:
                                newAttribute = attribute.split(';')
                                suffix,identity = newAttribute[5].split('ltr_identity=')
                                identity = identity.strip()
                                if identity != 'NA':
                                    identity = float(identity)
                                    dist = 1-identity
                                    time = (float(dist) / (2*lambda_value))
                                    time_in_thousands = float(time) / thousandsTimeDivisor
                                    time_in_millions = float(time) / millionsTimeDivisor
                                    if 'Gypsy' not in ty3TimeDictFullIntact:
                                        ty3TimeDictFullIntact['Gypsy'] = {}
                                    if chrID not in ty3TimeDictFullIntact['Gypsy']:
                                        ty3TimeDictFullIntact['Gypsy'][chrID] = {}
                                    if assemblyID not in ty3TimeDictFullIntact['Gypsy'][chrID]:
                                        ty3TimeDictFullIntact['Gypsy'][chrID][assemblyID] = []
                                    ty3TimeDictFullIntact['Gypsy'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Gypsy' not in ty3TimeDictSubsetIntact:
                                            ty3TimeDictSubsetIntact['Gypsy'] = {}
                                        if chrID not in ty3TimeDictSubsetIntact['Gypsy']:
                                            ty3TimeDictSubsetIntact['Gypsy'][chrID] = {}
                                        if assemblyID not in ty3TimeDictSubsetIntact['Gypsy'][chrID]:
                                            ty3TimeDictSubsetIntact['Gypsy'][chrID][assemblyID] = []
                                        ty3TimeDictSubsetIntact['Gypsy'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if 'LTR_retrotransposon' in feature and 'Copia' not in feature and 'Gypsy' not in feature:
                            if 'Method=homology' in attribute:
                                newAttribute = attribute.split(';')
                                suffix,identity = newAttribute[4].split('Identity=')
                                identity = identity.strip()
                                if identity != 'NA':
                                    identity = float(identity)
                                    dist = 1-identity
                                    time = (float(dist) / (2*lambda_value))
                                    time_in_thousands = float(time) / thousandsTimeDivisor
                                    time_in_millions = float(time) / millionsTimeDivisor
                                    if 'unknown' not in unknownTimeDictFullNonIntact:
                                        unknownTimeDictFullNonIntact['unknown'] = {}
                                    if chrID not in unknownTimeDictFullNonIntact['unknown']:
                                        unknownTimeDictFullNonIntact['unknown'][chrID] = {}
                                    if assemblyID not in unknownTimeDictFullNonIntact['unknown'][chrID]:
                                        unknownTimeDictFullNonIntact['unknown'][chrID][assemblyID] = []
                                    unknownTimeDictFullNonIntact['unknown'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'unknown' not in unknownTimeDictSubsetNonIntact:
                                            unknownTimeDictSubsetNonIntact['unknown'] = {}
                                        if chrID not in unknownTimeDictSubsetNonIntact['unknown']:
                                            unknownTimeDictSubsetNonIntact['unknown'][chrID] = {}
                                        if assemblyID not in unknownTimeDictSubsetNonIntact['unknown'][chrID]:
                                            unknownTimeDictSubsetNonIntact['unknown'][chrID][assemblyID] = []
                                        unknownTimeDictSubsetNonIntact['unknown'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                            if 'Method=structural' in attribute:
                                newAttribute = attribute.split(';')
                                suffix,identity = newAttribute[5].split('ltr_identity=')
                                identity = identity.strip()
                                if identity != 'NA':
                                    identity = float(identity)
                                    dist = 1-identity
                                    time = (float(dist) / (2*lambda_value))
                                    time_in_thousands = float(time) / thousandsTimeDivisor
                                    time_in_millions = float(time) / millionsTimeDivisor
                                    if 'unknown' not in unknownTimeDictFullIntact:
                                        unknownTimeDictFullIntact['unknown'] = {}
                                    if chrID not in unknownTimeDictFullIntact['unknown']:
                                        unknownTimeDictFullIntact['unknown'][chrID] = {}
                                    if assemblyID not in unknownTimeDictFullIntact['unknown'][chrID]:
                                        unknownTimeDictFullIntact['unknown'][chrID][assemblyID] = []
                                    unknownTimeDictFullIntact['unknown'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'unknown' not in unknownTimeDictSubsetIntact:
                                            unknownTimeDictSubsetIntact['unknown'] = {}
                                        if chrID not in unknownTimeDictSubsetIntact['unknown']:
                                            unknownTimeDictSubsetIntact['unknown'][chrID] = {}
                                        if assemblyID not in unknownTimeDictSubsetIntact['unknown'][chrID]:
                                            unknownTimeDictSubsetIntact['unknown'][chrID][assemblyID] = []
                                        unknownTimeDictSubsetIntact['unknown'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if feature == 'CACTA_TIR_transposon':
                            newAttribute = attribute.split(';')
                            suffix,identity = newAttribute[4].split('Identity=')
                            identity = identity.strip()
                            if identity != 'NA':
                                identity = float(identity)
                                dist = 1-identity
                                time = (float(dist) / (2*lambda_value))
                                time_in_thousands = float(time) / thousandsTimeDivisor
                                time_in_millions = float(time) / millionsTimeDivisor
                                if 'Method=homology' in attribute:
                                    if 'CACTA' not in CACTATimeDictFullNonIntact:
                                        CACTATimeDictFullNonIntact['CACTA'] = {}
                                    if chrID not in CACTATimeDictFullNonIntact['CACTA']:
                                        CACTATimeDictFullNonIntact['CACTA'][chrID] = {}
                                    if assemblyID not in CACTATimeDictFullNonIntact['CACTA'][chrID]:
                                        CACTATimeDictFullNonIntact['CACTA'][chrID][assemblyID] = []
                                    CACTATimeDictFullNonIntact['CACTA'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'CACTA' not in CACTATimeDictSubsetNonIntact:
                                            CACTATimeDictSubsetNonIntact['CACTA'] = {}
                                        if chrID not in CACTATimeDictSubsetNonIntact['CACTA']:
                                            CACTATimeDictSubsetNonIntact['CACTA'][chrID] = {}
                                        if assemblyID not in CACTATimeDictSubsetNonIntact['CACTA'][chrID]:
                                            CACTATimeDictSubsetNonIntact['CACTA'][chrID][assemblyID] = []
                                        CACTATimeDictSubsetNonIntact['CACTA'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                                if 'Method=structural' in attribute:
                                    if 'CACTA' not in CACTATimeDictFullIntact:
                                        CACTATimeDictFullIntact['CACTA'] = {}
                                    if chrID not in CACTATimeDictFullIntact['CACTA']:
                                        CACTATimeDictFullIntact['CACTA'][chrID] = {}
                                    if assemblyID not in CACTATimeDictFullIntact['CACTA'][chrID]:
                                        CACTATimeDictFullIntact['CACTA'][chrID][assemblyID] = []
                                    CACTATimeDictFullIntact['CACTA'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'CACTA' not in CACTATimeDictSubsetIntact:
                                            CACTATimeDictSubsetIntact['CACTA'] = {}
                                        if chrID not in CACTATimeDictSubsetIntact['CACTA']:
                                            CACTATimeDictSubsetIntact['CACTA'][chrID] = {}
                                        if assemblyID not in CACTATimeDictSubsetIntact['CACTA'][chrID]:
                                            CACTATimeDictSubsetIntact['CACTA'][chrID][assemblyID] = []
                                        CACTATimeDictSubsetIntact['CACTA'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if feature == 'Mutator_TIR_transposon':
                            newAttribute = attribute.split(';')
                            suffix,identity = newAttribute[4].split('Identity=')
                            identity = identity.strip()
                            if identity != 'NA':
                                identity = float(identity)
                                dist = 1-identity
                                time = (float(dist) / (2*lambda_value))
                                time_in_thousands = float(time) / thousandsTimeDivisor
                                time_in_millions = float(time) / millionsTimeDivisor
                                if 'Method=homology' in attribute:
                                    if 'Mutator' not in MutatorTimeDictFullNonIntact:
                                        MutatorTimeDictFullNonIntact['Mutator'] = {}
                                    if chrID not in MutatorTimeDictFullNonIntact['Mutator']:
                                        MutatorTimeDictFullNonIntact['Mutator'][chrID] = {}
                                    if assemblyID not in MutatorTimeDictFullNonIntact['Mutator'][chrID]:
                                        MutatorTimeDictFullNonIntact['Mutator'][chrID][assemblyID] = []
                                    MutatorTimeDictFullNonIntact['Mutator'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Mutator' not in MutatorTimeDictSubsetNonIntact:
                                            MutatorTimeDictSubsetNonIntact['Mutator'] = {}
                                        if chrID not in MutatorTimeDictSubsetNonIntact['Mutator']:
                                            MutatorTimeDictSubsetNonIntact['Mutator'][chrID] = {}
                                        if assemblyID not in MutatorTimeDictSubsetNonIntact['Mutator'][chrID]:
                                            MutatorTimeDictSubsetNonIntact['Mutator'][chrID][assemblyID] = []
                                        MutatorTimeDictSubsetNonIntact['Mutator'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                                if 'Method=structural' in attribute:
                                    if 'Mutator' not in MutatorTimeDictFullIntact:
                                        MutatorTimeDictFullIntact['Mutator'] = {}
                                    if chrID not in MutatorTimeDictFullIntact['Mutator']:
                                        MutatorTimeDictFullIntact['Mutator'][chrID] = {}
                                    if assemblyID not in MutatorTimeDictFullIntact['Mutator'][chrID]:
                                        MutatorTimeDictFullIntact['Mutator'][chrID][assemblyID] = []
                                    MutatorTimeDictFullIntact['Mutator'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Mutator' not in MutatorTimeDictSubsetIntact:
                                            MutatorTimeDictSubsetIntact['Mutator'] = {}
                                        if chrID not in MutatorTimeDictSubsetIntact['Mutator']:
                                            MutatorTimeDictSubsetIntact['Mutator'][chrID] = {}
                                        if assemblyID not in MutatorTimeDictSubsetIntact['Mutator'][chrID]:
                                            MutatorTimeDictSubsetIntact['Mutator'][chrID][assemblyID] = []
                                        MutatorTimeDictSubsetIntact['Mutator'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if feature == 'PIF_Harbinger_TIR_transposon':
                            newAttribute = attribute.split(';')
                            suffix,identity = newAttribute[4].split('Identity=')
                            identity = identity.strip()
                            if identity != 'NA':
                                identity = float(identity)
                                dist = 1-identity
                                time = (float(dist) / (2*lambda_value))
                                time_in_thousands = float(time) / thousandsTimeDivisor
                                time_in_millions = float(time) / millionsTimeDivisor
                                if 'Method=homology' in attribute:
                                    if 'Harbinger' not in HarbingerTimeDictFullNonIntact:
                                        HarbingerTimeDictFullNonIntact['Harbinger'] = {}
                                    if chrID not in HarbingerTimeDictFullNonIntact['Harbinger']:
                                        HarbingerTimeDictFullNonIntact['Harbinger'][chrID] = {}
                                    if assemblyID not in HarbingerTimeDictFullNonIntact['Harbinger'][chrID]:
                                        HarbingerTimeDictFullNonIntact['Harbinger'][chrID][assemblyID] = []
                                    HarbingerTimeDictFullNonIntact['Harbinger'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Harbinger' not in HarbingerTimeDictSubsetNonIntact:
                                            HarbingerTimeDictSubsetNonIntact['Harbinger'] = {}
                                        if chrID not in HarbingerTimeDictSubsetNonIntact['Harbinger']:
                                            HarbingerTimeDictSubsetNonIntact['Harbinger'][chrID] = {}
                                        if assemblyID not in HarbingerTimeDictSubsetNonIntact['Harbinger'][chrID]:
                                            HarbingerTimeDictSubsetNonIntact['Harbinger'][chrID][assemblyID] = []
                                        HarbingerTimeDictSubsetNonIntact['Harbinger'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                                if 'Method=structural' in attribute:
                                    if 'Harbinger' not in HarbingerTimeDictFullIntact:
                                        HarbingerTimeDictFullIntact['Harbinger'] = {}
                                    if chrID not in HarbingerTimeDictFullIntact['Harbinger']:
                                        HarbingerTimeDictFullIntact['Harbinger'][chrID] = {}
                                    if assemblyID not in HarbingerTimeDictFullIntact['Harbinger'][chrID]:
                                        HarbingerTimeDictFullIntact['Harbinger'][chrID][assemblyID] = []
                                    HarbingerTimeDictFullIntact['Harbinger'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Harbinger' not in HarbingerTimeDictSubsetIntact:
                                            HarbingerTimeDictSubsetIntact['Harbinger'] = {}
                                        if chrID not in HarbingerTimeDictSubsetIntact['Harbinger']:
                                            HarbingerTimeDictSubsetIntact['Harbinger'][chrID] = {}
                                        if assemblyID not in HarbingerTimeDictSubsetIntact['Harbinger'][chrID]:
                                            HarbingerTimeDictSubsetIntact['Harbinger'][chrID][assemblyID] = []
                                        HarbingerTimeDictSubsetIntact['Harbinger'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if feature == 'Tc1_Mariner_TIR_transposon':
                            newAttribute = attribute.split(';')
                            suffix,identity = newAttribute[4].split('Identity=')
                            identity = identity.strip()
                            if identity != 'NA':
                                identity = float(identity)
                                dist = 1-identity
                                time = (float(dist) / (2*lambda_value))
                                time_in_thousands = float(time) / thousandsTimeDivisor
                                time_in_millions = float(time) / millionsTimeDivisor
                                if 'Method=homology' in attribute:
                                    if 'Mariner' not in MarinerTimeDictFullNonIntact:
                                        MarinerTimeDictFullNonIntact['Mariner'] = {}
                                    if chrID not in MarinerTimeDictFullNonIntact['Mariner']:
                                        MarinerTimeDictFullNonIntact['Mariner'][chrID] = {}
                                    if assemblyID not in MarinerTimeDictFullNonIntact['Mariner'][chrID]:
                                        MarinerTimeDictFullNonIntact['Mariner'][chrID][assemblyID] = []
                                    MarinerTimeDictFullNonIntact['Mariner'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Mariner' not in MarinerTimeDictSubsetNonIntact:
                                            MarinerTimeDictSubsetNonIntact['Mariner'] = {}
                                        if chrID not in MarinerTimeDictSubsetNonIntact['Mariner']:
                                            MarinerTimeDictSubsetNonIntact['Mariner'][chrID] = {}
                                        if assemblyID not in MarinerTimeDictSubsetNonIntact['Mariner'][chrID]:
                                            MarinerTimeDictSubsetNonIntact['Mariner'][chrID][assemblyID] =     []
                                        MarinerTimeDictSubsetNonIntact['Mariner'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                                if 'Method=structural' in attribute:
                                    if 'Mariner' not in MarinerTimeDictFullIntact:
                                        MarinerTimeDictFullIntact['Mariner'] = {}
                                    if chrID not in MarinerTimeDictFullIntact['Mariner']:
                                        MarinerTimeDictFullIntact['Mariner'][chrID] = {}
                                    if assemblyID not in MarinerTimeDictFullIntact['Mariner'][chrID]:
                                        MarinerTimeDictFullIntact['Mariner'][chrID][assemblyID] = []
                                    MarinerTimeDictFullIntact['Mariner'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'Mariner' not in MarinerTimeDictSubsetIntact:
                                            MarinerTimeDictSubsetIntact['Mariner'] = {}
                                        if chrID not in MarinerTimeDictSubsetIntact['Mariner']:
                                            MarinerTimeDictSubsetIntact['Mariner'][chrID] = {}
                                        if assemblyID not in MarinerTimeDictSubsetIntact['Mariner'][chrID]:
                                            MarinerTimeDictSubsetIntact['Mariner'][chrID][assemblyID] = []
                                        MarinerTimeDictSubsetIntact['Mariner'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if feature == 'hAT_TIR_transposon':
                            newAttribute = attribute.split(';')
                            suffix,identity = newAttribute[4].split('Identity=')
                            identity = identity.strip()
                            if identity != 'NA':
                                identity = float(identity)
                                dist = 1-identity
                                time = (float(dist) / (2*lambda_value))
                                time_in_thousands = float(time) / thousandsTimeDivisor
                                time_in_millions = float(time) / millionsTimeDivisor
                                if 'Method=homology' in attribute:
                                    if 'hAT' not in hATTimeDictFullNonIntact:
                                        hATTimeDictFullNonIntact['hAT'] = {}
                                    if chrID not in hATTimeDictFullNonIntact['hAT']:
                                        hATTimeDictFullNonIntact['hAT'][chrID] = {}
                                    if assemblyID not in hATTimeDictFullNonIntact['hAT'][chrID]:
                                        hATTimeDictFullNonIntact['hAT'][chrID][assemblyID] = []
                                    hATTimeDictFullNonIntact['hAT'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'hAT' not in hATTimeDictSubsetNonIntact:
                                            hATTimeDictSubsetNonIntact['hAT'] = {}
                                        if chrID not in hATTimeDictSubsetNonIntact['hAT']:
                                            hATTimeDictSubsetNonIntact['hAT'][chrID] = {}
                                        if assemblyID not in hATTimeDictSubsetNonIntact['hAT'][chrID]:
                                            hATTimeDictSubsetNonIntact['hAT'][chrID][assemblyID] = []
                                        hATTimeDictSubsetNonIntact['hAT'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                                if 'Method=structural' in attribute:
                                    if 'hAT' not in hATTimeDictFullIntact:
                                        hATTimeDictFullIntact['hAT'] = {}
                                    if chrID not in hATTimeDictFullIntact['hAT']:
                                        hATTimeDictFullIntact['hAT'][chrID] = {}
                                    if assemblyID not in hATTimeDictFullIntact['hAT'][chrID]:
                                        hATTimeDictFullIntact['hAT'][chrID][assemblyID] = []
                                    hATTimeDictFullIntact['hAT'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'hAT' not in hATTimeDictSubsetIntact:
                                            hATTimeDictSubsetIntact['hAT'] = {}
                                        if chrID not in hATTimeDictSubsetIntact['hAT']:
                                            hATTimeDictSubsetIntact['hAT'][chrID] = {}
                                        if assemblyID not in hATTimeDictSubsetIntact['hAT'][chrID]:
                                            hATTimeDictSubsetIntact['hAT'][chrID][assemblyID] = []
                                        hATTimeDictSubsetIntact['hAT'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
                        if feature == 'helitron':
                            newAttribute = attribute.split(';')
                            suffix,identity = newAttribute[4].split('Identity=')
                            identity = identity.strip()
                            if identity != 'NA':
                                identity = float(identity)
                                dist = 1-identity
                                time = (float(dist) / (2*lambda_value))
                                time_in_thousands = float(time) / thousandsTimeDivisor
                                time_in_millions = float(time) / millionsTimeDivisor
                                if 'Method=homology' in attribute:
                                    if 'helitron' not in helitronTimeDictFullNonIntact:
                                        helitronTimeDictFullNonIntact['helitron'] = {}
                                    if chrID not in helitronTimeDictFullNonIntact['helitron']:
                                        helitronTimeDictFullNonIntact['helitron'][chrID] = {}
                                    if assemblyID not in helitronTimeDictFullNonIntact['helitron'][chrID]:
                                        helitronTimeDictFullNonIntact['helitron'][chrID][assemblyID] = []
                                    helitronTimeDictFullNonIntact['helitron'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'helitron' not in helitronTimeDictSubsetNonIntact:
                                            helitronTimeDictSubsetNonIntact['helitron'] = {}
                                        if chrID not in helitronTimeDictSubsetNonIntact['helitron']:
                                            helitronTimeDictSubsetNonIntact['helitron'][chrID] = {}
                                        if assemblyID not in helitronTimeDictSubsetNonIntact['helitron'][chrID]:
                                            helitronTimeDictSubsetNonIntact['helitron'][chrID][assemblyID] = []
                                        helitronTimeDictSubsetNonIntact['helitron'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("NonIntact",feature,chrID,assemblyID,time_in_thousands))
                                if 'Method=structural' in attribute:
                                    if 'helitron' not in helitronTimeDictFullIntact:
                                        helitronTimeDictFullIntact['helitron'] = {}
                                    if chrID not in helitronTimeDictFullIntact['helitron']:
                                        helitronTimeDictFullIntact['helitron'][chrID] = {}
                                    if assemblyID not in helitronTimeDictFullIntact['helitron'][chrID]:
                                        helitronTimeDictFullIntact['helitron'][chrID][assemblyID] = []
                                    helitronTimeDictFullIntact['helitron'][chrID][assemblyID].append(time_in_millions)
                                    OUT_full.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_millions))
                                    if identity >= distanceThreshold:
                                        if 'helitron' not in helitronTimeDictSubsetIntact:
                                            helitronTimeDictSubsetIntact['helitron'] = {}
                                        if chrID not in helitronTimeDictSubsetIntact['helitron']:
                                            helitronTimeDictSubsetIntact['helitron'][chrID] = {}
                                        if assemblyID not in helitronTimeDictSubsetIntact['helitron'][chrID]:
                                            helitronTimeDictSubsetIntact['helitron'][chrID][assemblyID] = []
                                        helitronTimeDictSubsetIntact['helitron'][chrID][assemblyID].append(time_in_thousands)
                                        OUT_subset.write("%s\t%s\t%s\t%s\t%s\n" % ("Intact",feature,chrID,assemblyID,time_in_thousands))
    return(CACTATimeDictFullNonIntact,MutatorTimeDictFullNonIntact,HarbingerTimeDictFullNonIntact,MarinerTimeDictFullNonIntact,hATTimeDictFullNonIntact,helitronTimeDictFullNonIntact,ty1TimeDictFullNonIntact,ty3TimeDictFullNonIntact,unknownTimeDictFullNonIntact,CACTATimeDictSubsetNonIntact,MutatorTimeDictSubsetNonIntact,HarbingerTimeDictSubsetNonIntact,MarinerTimeDictSubsetNonIntact,hATTimeDictSubsetNonIntact,helitronTimeDictSubsetNonIntact,ty1TimeDictSubsetNonIntact,ty3TimeDictSubsetNonIntact,unknownTimeDictSubsetNonIntact,CACTATimeDictFullIntact,MutatorTimeDictFullIntact,HarbingerTimeDictFullIntact,MarinerTimeDictFullIntact,hATTimeDictFullIntact,helitronTimeDictFullIntact,ty1TimeDictFullIntact,ty3TimeDictFullIntact,unknownTimeDictFullIntact, ty1TimeDictSubsetIntact, ty3TimeDictSubsetIntact, unknownTimeDictSubsetIntact, CACTATimeDictSubsetIntact, MutatorTimeDictSubsetIntact, HarbingerTimeDictSubsetIntact, MarinerTimeDictSubsetIntact, hATTimeDictSubsetIntact, helitronTimeDictSubsetIntact)


def createTEDataLists(dataDict):
    fullList = []
    for teType in dataDict:
        for chrID in dataDict[teType]:
            # initialize figures here
            for assemblyID in dataDict[teType][chrID]:
                for time in dataDict[teType][chrID][assemblyID]:
                    fullList.append(time)
    return(fullList)


# get the data into chr-level data structure
def createChrLevelDataDict(dataDict):
    chrLevelDataDict = {}
    for ltrType in dataDict:
        if ltrType not in chrLevelDataDict:
            chrLevelDataDict[ltrType] = {}
        for chrID in dataDict[ltrType]:
            if chrID not in chrLevelDataDict[ltrType]:
                chrLevelDataDict[ltrType][chrID] = []
            for assemblyID in dataDict[ltrType][chrID]:
                for time in dataDict[ltrType][chrID][assemblyID]:
                    chrLevelDataDict[ltrType][chrID].append(time)
    return(chrLevelDataDict)


def createCombinedDataHist(NonIntactCACTAList,NonIntactMutatorList,NonIntactHarbingerList,NonIntactMarinerList,NonIntact_hATList,NonIntact_helitronList,NonIntactTy1List,NonIntactTy3List,NonIntactUnknownList,IntactCACTAList,IntactMutatorList,IntactHarbingerList,IntactMarinerList,Intact_hATList,Intact_helitronList,IntactTy1List,IntactTy3List,IntactUnknownList,timeLabel,figLabel,binWidth,teSet):
    #plt.rcParams["figure.figsize"] = (5,4)
    colorMap = {'CACTA':'#000000','Mutator':'#E69F00','PIF_Harbinger':'#0072B2','Tc1_Mariner':'#009E73','hAT_TIR':'#F0E442','helitron':'#CC79A7','Ty1':'#542788', 'unknown':'#998ec3', 'Ty3':'#b35806'}
    # 'Ty3':'#f1a340'
    # binWidth = 0.5
    binWidth = binWidth

    completeNonIntactDataList = NonIntactCACTAList+NonIntactMutatorList+NonIntactHarbingerList+NonIntactMarinerList+NonIntact_hATList+NonIntact_helitronList+NonIntactTy1List+NonIntactTy3List+NonIntactUnknownList
    ### fig, ax1 = plt.subplots(figsize=(5,4))
    # fig, ax1 = plt.subplots(figsize=(7.5,6))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[10,6])
    bins = np.arange(min(completeNonIntactDataList),max(completeNonIntactDataList),binWidth)
    counts1, bins1, bars1 = ax1.hist(NonIntactTy3List, bins = bins, histtype = 'step', color = colorMap['Ty3'], label='Ty3-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts2, bins2, bars2 = ax1.hist(NonIntactTy1List, bins = bins, histtype = 'step', color = colorMap['Ty1'], label='Ty1-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts3, bins3, bars3 = ax1.hist(NonIntactUnknownList, bins = bins, histtype = 'step', color = colorMap['unknown'], label='Unknown-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')

    counts4, bins4, bars4 = ax1.hist(NonIntactMutatorList, bins = bins, histtype = 'step', color = colorMap['Mutator'], label='Mutator', linewidth=2, alpha=0.75)
    counts5, bins5, bars5 = ax1.hist(NonIntactHarbingerList, bins = bins, histtype = 'step', color = colorMap['PIF_Harbinger'], label='Harbinger', linewidth=2, alpha=0.75)
    counts6, bins6, bars6 = ax1.hist(NonIntactMarinerList, bins = bins, histtype = 'step', color = colorMap['Tc1_Mariner'], label='Mariner', linewidth=2, alpha=0.75)
    counts7, bins7, bars7 = ax1.hist(NonIntact_hATList, bins = bins, histtype = 'step', color = colorMap['hAT_TIR'], label='hAT', linewidth=2, alpha=0.75)
    counts8, bins8, bars8 = ax1.hist(NonIntact_helitronList, bins = bins, histtype = 'step', color = colorMap['helitron'], label='Helitron', linewidth=2, alpha=0.75)
    counts9, bins9, bars9 = ax1.hist(NonIntactCACTAList, bins = bins, histtype = 'step', color = colorMap['CACTA'], label='CACTA', linewidth=2, alpha=0.75)
    ax1.margins(0)
    ax1.set_ylabel('Transposable element count (all scaffolded assemblies)')

    #left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
    #ax2 = fig.add_axes([left, bottom, width, height])
    #axins2 = inset_axes(ax2, width="30%", height="50%")
    completeIntactDataList = IntactCACTAList+IntactMutatorList+IntactHarbingerList+IntactMarinerList+Intact_hATList+Intact_helitronList+IntactTy1List+IntactTy3List+IntactUnknownList
    bins = np.arange(min(completeIntactDataList),max(completeIntactDataList),binWidth)
    counts1, bins1, bars1 = ax2.hist(IntactTy3List, bins = bins, histtype = 'step', color = colorMap['Ty3'], label='Ty3-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts2, bins2, bars2 = ax2.hist(IntactTy1List, bins = bins, histtype = 'step', color = colorMap['Ty1'], label='Ty1-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts3, bins3, bars3 = ax2.hist(IntactUnknownList, bins = bins, histtype = 'step', color = colorMap['unknown'], label='Unknown-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')

    counts4, bins4, bars4 = ax2.hist(IntactMutatorList, bins = bins, histtype = 'step', color = colorMap['Mutator'], label='Mutator', linewidth=2, alpha=0.75)
    counts5, bins5, bars5 = ax2.hist(IntactHarbingerList, bins = bins, histtype = 'step', color = colorMap['PIF_Harbinger'], label='Harbinger', linewidth=2, alpha=0.75)
    counts6, bins6, bars6 = ax2.hist(IntactMarinerList, bins = bins, histtype = 'step', color = colorMap['Tc1_Mariner'], label='Mariner', linewidth=2, alpha=0.75)
    counts7, bins7, bars7 = ax2.hist(Intact_hATList, bins = bins, histtype = 'step', color = colorMap['hAT_TIR'], label='hAT', linewidth=2, alpha=0.75)
    counts8, bins8, bars8 = ax2.hist(Intact_helitronList, bins = bins, histtype = 'step', color = colorMap['helitron'], label='Helitron', linewidth=2, alpha=0.75)
    counts9, bins9, bars9 = ax2.hist(IntactCACTAList, bins = bins, histtype = 'step', color = colorMap['CACTA'], label='CACTA', linewidth=2, alpha=0.75)
    ax2.margins(0)

    ax1.set_title('Non-intact TEs')
    ax2.set_title('Intact TEs')
    fig.supxlabel('Time (' + timeLabel + ')')
    #plt.xlabel('Time (' + timeLabel + ')')
    #ax1.set_xlabel('Time (' + timeLabel + ')')
    #plt.margins(0)
    legendLabelTy1 = mlines.Line2D([], [], color=colorMap['Ty1'], label='Ty1-LTR', linestyle='dashdot', linewidth=2)
    legendLabelTy3 = mlines.Line2D([], [], color=colorMap['Ty3'], label='Ty3-LTR', linestyle='dashdot', linewidth=2)
    legendLabelUnknown = mlines.Line2D([], [], color=colorMap['unknown'], label='Unknown-LTR', linestyle='dashdot', linewidth=2)
    #legendLabelCACTA = mlines.Line2D([], [], color=colorMap['CACTA'], marker='s',linestyle="None",markersize=10, label='CACTA')
    legendLabelCACTA = mlines.Line2D([], [], color=colorMap['CACTA'], label='CACTA', linestyle='solid', linewidth=2)
    legendLabelMutator = mlines.Line2D([], [], color=colorMap['Mutator'], label='Mutator', linestyle='solid', linewidth=2)
    legendLabelHarbinger = mlines.Line2D([], [], color=colorMap['PIF_Harbinger'], label='Harbinger', linestyle='solid', linewidth=2)
    legendLabelMariner = mlines.Line2D([], [], color=colorMap['Tc1_Mariner'], label='Mariner', linestyle='solid', linewidth=2)
    legendLabel_hAT = mlines.Line2D([], [], color=colorMap['hAT_TIR'], label='hAT', linestyle='solid', linewidth=2)
    legendLabelHelitron = mlines.Line2D([], [], color=colorMap['helitron'], label='Helitron', linestyle='solid', linewidth=2)
    
    ax2.legend(handles=[legendLabelTy3,legendLabelTy1,legendLabelUnknown,legendLabelCACTA,
                        legendLabelMutator,legendLabelHarbinger,legendLabelMariner,
                        legendLabel_hAT,legendLabelHelitron], frameon=False, ncol=2, loc='upper right')
    plt.savefig(figLabel + '_allChroms_' + teSet + 'TEDistances.svg')
    plt.savefig(figLabel + '_allChroms_' + teSet + 'TEDistances.png', dpi=600)
    plt.close()
    

def createCombinedDataHistWithInset(fullNonIntactCACTAList,fullNonIntactMutatorList,fullNonIntactHarbingerList,fullNonIntactMarinerList,fullNonIntact_hATList,fullNonIntact_helitronList,fullNonIntactTy1List,fullNonIntactTy3List,fullNonIntactUnknownList,fullIntactCACTAList,fullIntactMutatorList,fullIntactHarbingerList,fullIntactMarinerList,fullIntact_hATList,fullIntact_helitronList,fullIntactTy1List,fullIntactTy3List,fullIntactUnknownList,full_timeLabel,full_figLabel,full_binWidth,full_teSet,subsetNonIntactCACTAList,subsetNonIntactMutatorList,subsetNonIntactHarbingerList,subsetNonIntactMarinerList,subsetNonIntact_hATList,subsetNonIntact_helitronList,subsetNonIntactTy1List,subsetNonIntactTy3List,subsetNonIntactUnknownList,subsetIntactCACTAList,subsetIntactMutatorList,subsetIntactHarbingerList,subsetIntactMarinerList,subsetIntact_hATList,subsetIntact_helitronList,subsetIntactTy1List,subsetIntactTy3List,subsetIntactUnknownList,subset_timeLabel,subset_figLabel,subset_binWidth,subset_teSet):
    colorMap = {'CACTA':'#000000','Mutator':'#E69F00','PIF_Harbinger':'#0072B2','Tc1_Mariner':'#009E73','hAT_TIR':'#F0E442','helitron':'#CC79A7','Ty1':'#542788', 'unknown':'#998ec3', 'Ty3':'#b35806'}
    
    binWidth = full_binWidth
    full_completeNonIntactDataList = fullNonIntactCACTAList+fullNonIntactMutatorList+fullNonIntactHarbingerList+fullNonIntactMarinerList+fullNonIntact_hATList+fullNonIntact_helitronList+fullNonIntactTy1List+fullNonIntactTy3List+fullNonIntactUnknownList
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[12,7])
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[10,6])
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[20,12])
    bins = np.arange(min(full_completeNonIntactDataList),max(full_completeNonIntactDataList),binWidth)
    counts1, bins1, bars1 = ax1.hist(fullNonIntactTy3List, bins = bins, histtype = 'step', color = colorMap['Ty3'], label='Ty3-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts2, bins2, bars2 = ax1.hist(fullNonIntactTy1List, bins = bins, histtype = 'step', color = colorMap['Ty1'], label='Ty1-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts3, bins3, bars3 = ax1.hist(fullNonIntactUnknownList, bins = bins, histtype = 'step', color = colorMap['unknown'], label='Unknown-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts4, bins4, bars4 = ax1.hist(fullNonIntactMutatorList, bins = bins, histtype = 'step', color = colorMap['Mutator'], label='Mutator', linewidth=2, alpha=0.75)
    counts5, bins5, bars5 = ax1.hist(fullNonIntactHarbingerList, bins = bins, histtype = 'step', color = colorMap['PIF_Harbinger'], label='Harbinger', linewidth=2, alpha=0.75)
    counts6, bins6, bars6 = ax1.hist(fullNonIntactMarinerList, bins = bins, histtype = 'step', color = colorMap['Tc1_Mariner'], label='Mariner', linewidth=2, alpha=0.75)
    counts7, bins7, bars7 = ax1.hist(fullNonIntact_hATList, bins = bins, histtype = 'step', color = colorMap['hAT_TIR'], label='hAT', linewidth=2, alpha=0.75)
    counts8, bins8, bars8 = ax1.hist(fullNonIntact_helitronList, bins = bins, histtype = 'step', color = colorMap['helitron'], label='Helitron', linewidth=2, alpha=0.75)
    counts9, bins9, bars9 = ax1.hist(fullNonIntactCACTAList, bins = bins, histtype = 'step', color = colorMap['CACTA'], label='CACTA', linewidth=2, alpha=0.75)
    ax1.margins(0)

    ax3 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax1, [0.73,0.7,0.25,0.25])
    # inset_axes
    
    # ip = inset_axes(ax1, [0.73,0.7,0.25,0.25])
    # InsetPosition(ax1, [% of x-axis coords, % of y-axis coords, % of width, % of height])
    # The following bounds the inset axes to a box with 20% of the parent axes's height and 40% of the width.
    # ip = InsetPosition(ax, [0.5, 0.1, 0.4, 0.2])
    ax3.set_axes_locator(ip)
    # mark_inset(ax1, ax3, loc1=2, loc2=4, fc="none", ec='0.5')
    ###
    binWidth = subset_binWidth
    subset_completeNonIntactDataList = subsetNonIntactCACTAList+subsetNonIntactMutatorList+subsetNonIntactHarbingerList+subsetNonIntactMarinerList+subsetNonIntact_hATList+subsetNonIntact_helitronList+subsetNonIntactTy1List+subsetNonIntactTy3List+subsetNonIntactUnknownList
    bins = np.arange(min(subset_completeNonIntactDataList),max(subset_completeNonIntactDataList),binWidth)
    counts1, bins1, bars1 = ax3.hist(subsetNonIntactTy3List, bins = bins, histtype = 'step', color = colorMap['Ty3'], label='Ty3-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts2, bins2, bars2 = ax3.hist(subsetNonIntactTy1List, bins = bins, histtype = 'step', color = colorMap['Ty1'], label='Ty1-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts3, bins3, bars3 = ax3.hist(subsetNonIntactUnknownList, bins = bins, histtype = 'step', color = colorMap['unknown'], label='Unknown-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts4, bins4, bars4 = ax3.hist(subsetNonIntactMutatorList, bins = bins, histtype = 'step', color = colorMap['Mutator'], label='Mutator', linewidth=2, alpha=0.75)
    counts5, bins5, bars5 = ax3.hist(subsetNonIntactHarbingerList, bins = bins, histtype = 'step', color = colorMap['PIF_Harbinger'], label='Harbinger', linewidth=2, alpha=0.75)
    counts6, bins6, bars6 = ax3.hist(subsetNonIntactMarinerList, bins = bins, histtype = 'step', color = colorMap['Tc1_Mariner'], label='Mariner', linewidth=2, alpha=0.75)
    counts7, bins7, bars7 = ax3.hist(subsetNonIntact_hATList, bins = bins, histtype = 'step', color = colorMap['hAT_TIR'], label='hAT', linewidth=2, alpha=0.75)
    counts8, bins8, bars8 = ax3.hist(subsetNonIntact_helitronList, bins = bins, histtype = 'step', color = colorMap['helitron'], label='Helitron', linewidth=2, alpha=0.75)
    counts9, bins9, bars9 = ax3.hist(subsetNonIntactCACTAList, bins = bins, histtype = 'step', color = colorMap['CACTA'], label='CACTA', linewidth=2, alpha=0.75)
    ax3.margins(0)
    ax3.set_xlabel('Time (' + subset_timeLabel + ')')
    ###

    full_completeIntactDataList = fullIntactCACTAList+fullIntactMutatorList+fullIntactHarbingerList+fullIntactMarinerList+fullIntact_hATList+fullIntact_helitronList+fullIntactTy1List+fullIntactTy3List+fullIntactUnknownList
    binWidth = full_binWidth
    bins = np.arange(min(full_completeIntactDataList),max(full_completeIntactDataList),binWidth)
    counts1, bins1, bars1 = ax2.hist(fullIntactTy3List, bins = bins, histtype = 'step', color = colorMap['Ty3'], label='Ty3-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts2, bins2, bars2 = ax2.hist(fullIntactTy1List, bins = bins, histtype = 'step', color = colorMap['Ty1'], label='Ty1-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts3, bins3, bars3 = ax2.hist(fullIntactUnknownList, bins = bins, histtype = 'step', color = colorMap['unknown'], label='Unknown-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts4, bins4, bars4 = ax2.hist(fullIntactMutatorList, bins = bins, histtype = 'step', color = colorMap['Mutator'], label='Mutator', linewidth=2, alpha=0.75)
    counts5, bins5, bars5 = ax2.hist(fullIntactHarbingerList, bins = bins, histtype = 'step', color = colorMap['PIF_Harbinger'], label='Harbinger', linewidth=2, alpha=0.75)
    counts6, bins6, bars6 = ax2.hist(fullIntactMarinerList, bins = bins, histtype = 'step', color = colorMap['Tc1_Mariner'], label='Mariner', linewidth=2, alpha=0.75)
    counts7, bins7, bars7 = ax2.hist(fullIntact_hATList, bins = bins, histtype = 'step', color = colorMap['hAT_TIR'], label='hAT', linewidth=2, alpha=0.75)
    counts8, bins8, bars8 = ax2.hist(fullIntact_helitronList, bins = bins, histtype = 'step', color = colorMap['helitron'], label='Helitron', linewidth=2, alpha=0.75)
    counts9, bins9, bars9 = ax2.hist(fullIntactCACTAList, bins = bins, histtype = 'step', color = colorMap['CACTA'], label='CACTA', linewidth=2, alpha=0.75)
    ax2.margins(0)

    ax4 = plt.axes([0,0,1,1])
    # Axes.inset_axes
    ip = InsetPosition(ax2, [0.73,0.7,0.25,0.25])
    # ip = inset_axes(ax2, [0.73,0.7,0.25,0.25])
    # ip = ax2.inset_axes([0.73,0.7,0.25,0.25])

    ax4.set_axes_locator(ip)
    
    subset_completeIntactDataList = subsetIntactCACTAList+subsetIntactMutatorList+subsetIntactHarbingerList+subsetIntactMarinerList+subsetIntact_hATList+subsetIntact_helitronList+subsetIntactTy1List+subsetIntactTy3List+subsetIntactUnknownList
    binWidth = subset_binWidth
    bins = np.arange(min(subset_completeIntactDataList),max(subset_completeIntactDataList),binWidth)
    counts1, bins1, bars1 = ax4.hist(subsetIntactTy3List, bins = bins, histtype = 'step', color = colorMap['Ty3'], label='Ty3-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts2, bins2, bars2 = ax4.hist(subsetIntactTy1List, bins = bins, histtype = 'step', color = colorMap['Ty1'], label='Ty1-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts3, bins3, bars3 = ax4.hist(subsetIntactUnknownList, bins = bins, histtype = 'step', color = colorMap['unknown'], label='Unknown-LTR', linewidth=2, alpha=0.75, linestyle='dashdot')
    counts4, bins4, bars4 = ax4.hist(subsetIntactMutatorList, bins = bins, histtype = 'step', color = colorMap['Mutator'], label='Mutator', linewidth=2, alpha=0.75)
    counts5, bins5, bars5 = ax4.hist(subsetIntactHarbingerList, bins = bins, histtype = 'step', color = colorMap['PIF_Harbinger'], label='Harbinger', linewidth=2, alpha=0.75)
    counts6, bins6, bars6 = ax4.hist(subsetIntactMarinerList, bins = bins, histtype = 'step', color = colorMap['Tc1_Mariner'], label='Mariner', linewidth=2, alpha=0.75)
    counts7, bins7, bars7 = ax4.hist(subsetIntact_hATList, bins = bins, histtype = 'step', color = colorMap['hAT_TIR'], label='hAT', linewidth=2, alpha=0.75)
    counts8, bins8, bars8 = ax4.hist(subsetIntact_helitronList, bins = bins, histtype = 'step', color = colorMap['helitron'], label='Helitron', linewidth=2, alpha=0.75)
    counts9, bins9, bars9 = ax4.hist(subsetIntactCACTAList, bins = bins, histtype = 'step', color = colorMap['CACTA'], label='CACTA', linewidth=2, alpha=0.75)
    ax4.margins(0)
    ax4.set_xlabel('Time (' + subset_timeLabel + ')')
    
    ax1.set_title('Non-intact TEs')
    ax2.set_title('Intact TEs')
    fig.supxlabel('Time (' + full_timeLabel + ')')
    legendLabelTy1 = mlines.Line2D([], [], color=colorMap['Ty1'], label='Ty1-LTR', linestyle='dashdot', linewidth=2)
    legendLabelTy3 = mlines.Line2D([], [], color=colorMap['Ty3'], label='Ty3-LTR', linestyle='dashdot', linewidth=2)
    legendLabelUnknown = mlines.Line2D([], [], color=colorMap['unknown'], label='Unknown-LTR', linestyle='dashdot', linewidth=2)
    legendLabelCACTA = mlines.Line2D([], [], color=colorMap['CACTA'], label='CACTA', linestyle='solid', linewidth=2)
    legendLabelMutator = mlines.Line2D([], [], color=colorMap['Mutator'], label='Mutator', linestyle='solid', linewidth=2)
    legendLabelHarbinger = mlines.Line2D([], [], color=colorMap['PIF_Harbinger'], label='Harbinger', linestyle='solid', linewidth=2)
    legendLabelMariner = mlines.Line2D([], [], color=colorMap['Tc1_Mariner'], label='Mariner', linestyle='solid', linewidth=2)
    legendLabel_hAT = mlines.Line2D([], [], color=colorMap['hAT_TIR'], label='hAT', linestyle='solid', linewidth=2)
    legendLabelHelitron = mlines.Line2D([], [], color=colorMap['helitron'], label='Helitron', linestyle='solid', linewidth=2)

    ax2.legend(handles=[legendLabelTy3,legendLabelTy1,legendLabelUnknown,legendLabelCACTA,
                        legendLabelMutator,legendLabelHarbinger,legendLabelMariner,
                        legendLabel_hAT,legendLabelHelitron], frameon=False, ncol=1, loc='center')
    plt.savefig(full_figLabel + '_allChroms_' + full_teSet + 'TEDistances.svg')
    plt.savefig(full_figLabel + '_allChroms_' + full_teSet + 'TEDistances.png', dpi=600)
    plt.close()
    

    
usage = "Usage: " + sys.argv[0] + " <te anno gff file list>\n"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

teAnnoGFFFileList = sys.argv[1]

teAnnoGFFDataList = readFileList(teAnnoGFFFileList)

CACTATimeDictFullNonIntact,MutatorTimeDictFullNonIntact,HarbingerTimeDictFullNonIntact,MarinerTimeDictFullNonIntact,hATTimeDictFullNonIntact,helitronTimeDictFullNonIntact,ty1TimeDictFullNonIntact,ty3TimeDictFullNonIntact,unknownTimeDictFullNonIntact,CACTATimeDictSubsetNonIntact,MutatorTimeDictSubsetNonIntact,HarbingerTimeDictSubsetNonIntact,MarinerTimeDictSubsetNonIntact,hATTimeDictSubsetNonIntact,helitronTimeDictSubsetNonIntact,ty1TimeDictSubsetNonIntact,ty3TimeDictSubsetNonIntact,unknownTimeDictSubsetNonIntact,CACTATimeDictFullIntact,MutatorTimeDictFullIntact,HarbingerTimeDictFullIntact,MarinerTimeDictFullIntact,hATTimeDictFullIntact,helitronTimeDictFullIntact,ty1TimeDictFullIntact,ty3TimeDictFullIntact,unknownTimeDictFullIntact,ty1TimeDictSubsetIntact, ty3TimeDictSubsetIntact, unknownTimeDictSubsetIntact, CACTATimeDictSubsetIntact, MutatorTimeDictSubsetIntact, HarbingerTimeDictSubsetIntact, MarinerTimeDictSubsetIntact, hATTimeDictSubsetIntact, helitronTimeDictSubsetIntact = readTEAnnoGFF(teAnnoGFFDataList)

# full data
fullNonIntactCACTAList = createTEDataLists(CACTATimeDictFullNonIntact)
fullNonIntactMutatorList = createTEDataLists(MutatorTimeDictFullNonIntact)
fullNonIntactHarbingerList = createTEDataLists(HarbingerTimeDictFullNonIntact)
fullNonIntactMarinerList = createTEDataLists(MarinerTimeDictFullNonIntact)
fullNonIntact_hATList = createTEDataLists(hATTimeDictFullNonIntact)
fullNonIntact_helitronList = createTEDataLists(helitronTimeDictFullNonIntact)
fullNonIntactTy1List = createTEDataLists(ty1TimeDictFullNonIntact)
fullNonIntactTy3List = createTEDataLists(ty3TimeDictFullNonIntact)
fullNonIntactUnknownList = createTEDataLists(unknownTimeDictFullNonIntact)

fullIntactCACTAList = createTEDataLists(CACTATimeDictFullIntact)
fullIntactMutatorList = createTEDataLists(MutatorTimeDictFullIntact)
fullIntactHarbingerList = createTEDataLists(HarbingerTimeDictFullIntact)
fullIntactMarinerList = createTEDataLists(MarinerTimeDictFullIntact)
fullIntact_hATList = createTEDataLists(hATTimeDictFullIntact)
fullIntact_helitronList = createTEDataLists(helitronTimeDictFullIntact)
fullIntactTy1List = createTEDataLists(ty1TimeDictFullIntact)
fullIntactTy3List = createTEDataLists(ty3TimeDictFullIntact)
fullIntactUnknownList = createTEDataLists(unknownTimeDictFullIntact)

# subset
subsetNonIntactCACTAList = createTEDataLists(CACTATimeDictSubsetNonIntact)
subsetNonIntactMutatorList = createTEDataLists(MutatorTimeDictSubsetNonIntact)
subsetNonIntactHarbingerList = createTEDataLists(HarbingerTimeDictSubsetNonIntact)
subsetNonIntactMarinerList = createTEDataLists(MarinerTimeDictSubsetNonIntact)
subsetNonIntact_hATList = createTEDataLists(hATTimeDictSubsetNonIntact)
subsetNonIntact_helitronList = createTEDataLists(helitronTimeDictSubsetNonIntact)
subsetNonIntactTy1List = createTEDataLists(ty1TimeDictSubsetNonIntact)
subsetNonIntactTy3List = createTEDataLists(ty3TimeDictSubsetNonIntact)
subsetNonIntactUnknownList = createTEDataLists(unknownTimeDictSubsetNonIntact)

subsetIntactCACTAList = createTEDataLists(CACTATimeDictSubsetIntact)
subsetIntactMutatorList = createTEDataLists(MutatorTimeDictSubsetIntact)
subsetIntactHarbingerList = createTEDataLists(HarbingerTimeDictSubsetIntact)
subsetIntactMarinerList = createTEDataLists(MarinerTimeDictSubsetIntact)
subsetIntact_hATList = createTEDataLists(hATTimeDictSubsetIntact)
subsetIntact_helitronList = createTEDataLists(helitronTimeDictSubsetIntact)
subsetIntactTy1List = createTEDataLists(ty1TimeDictSubsetIntact)
subsetIntactTy3List = createTEDataLists(ty3TimeDictSubsetIntact)
subsetIntactUnknownList = createTEDataLists(unknownTimeDictSubsetIntact)

#createCombinedDataHist(fullNonIntactCACTAList,fullNonIntactMutatorList,fullNonIntactHarbingerList,fullNonIntactMarinerList,fullNonIntact_hATList,fullNonIntact_helitronList,fullNonIntactTy1List,fullNonIntactTy3List,fullNonIntactUnknownList,fullIntactCACTAList,fullIntactMutatorList,fullIntactHarbingerList,fullIntactMarinerList,fullIntact_hATList,fullIntact_helitronList,fullIntactTy1List,fullIntactTy3List,fullIntactUnknownList,'million years ago','allTEs',0.5,'all')

#createCombinedDataHist(subsetNonIntactCACTAList,subsetNonIntactMutatorList,subsetNonIntactHarbingerList,subsetNonIntactMarinerList,subsetNonIntact_hATList,subsetNonIntact_helitronList,subsetNonIntactTy1List,subsetNonIntactTy3List,subsetNonIntactUnknownList,subsetIntactCACTAList,subsetIntactMutatorList,subsetIntactHarbingerList,subsetIntactMarinerList,subsetIntact_hATList,subsetIntact_helitronList,subsetIntactTy1List,subsetIntactTy3List,subsetIntactUnknownList,'thousand years ago','allTEs',10,'subset')

createCombinedDataHistWithInset(fullNonIntactCACTAList,fullNonIntactMutatorList,fullNonIntactHarbingerList,fullNonIntactMarinerList,fullNonIntact_hATList,fullNonIntact_helitronList,fullNonIntactTy1List,fullNonIntactTy3List,fullNonIntactUnknownList,fullIntactCACTAList,fullIntactMutatorList,fullIntactHarbingerList,fullIntactMarinerList,fullIntact_hATList,fullIntact_helitronList,fullIntactTy1List,fullIntactTy3List,fullIntactUnknownList,'mya','allTEs',0.5,'all',subsetNonIntactCACTAList,subsetNonIntactMutatorList,subsetNonIntactHarbingerList,subsetNonIntactMarinerList,subsetNonIntact_hATList,subsetNonIntact_helitronList,subsetNonIntactTy1List,subsetNonIntactTy3List,subsetNonIntactUnknownList,subsetIntactCACTAList,subsetIntactMutatorList,subsetIntactHarbingerList,subsetIntactMarinerList,subsetIntact_hATList,subsetIntact_helitronList,subsetIntactTy1List,subsetIntactTy3List,subsetIntactUnknownList,'kya','allTEs',10,'subset')

# million years ago
# thousand years ago

#fullTy1List = fullNonIntactTy1List + fullCopiaList
#fullTy3List = fullNonIntactTy3List + fullGypsyList
#fullUnknownList = fullNonIntactUnknownList + fullOtherList

'''
createCombinedFullDataHist(fullTy1List,fullTy3List,fullUnknownList,
                           fullCACTAList,fullMutatorList,fullHarbingerList,fullMarinerList,full_hATList,full_helitronList,
                           'million years ago','allTEs',0.5)
'''

'''
createCombinedFullDataHist(fullGypsyList,fullCopiaList,fullOtherList,
                   fullCACTAList,fullMutatorList,fullHarbingerList,fullMarinerList,full_hATList,full_helitronList,
                           'million years ago','allTEs',0.5)
'''

