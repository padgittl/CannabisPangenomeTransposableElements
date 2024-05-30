#!/bin/python
import sys, os, re
import numpy as np


###############
# SUBROUTINES #
###############


def readBedFile(bedFile,specificAssemblyID,moleculeID):
    assemblyIDs = {}
    with open(bedFile,'r') as F:
        for line in F:
            scaffoldID,start,stop,synthaseModelID,score,strand = line.strip().split('\t')
            assemblyID,chromID = scaffoldID.split('.')
            if assemblyID not in assemblyIDs:
                assemblyIDs[assemblyID] = 1
    if specificAssemblyID in assemblyIDs:
        OUT = open(specificAssemblyID + '_' + moleculeID + '_fullhits.bed','w')
        with open(bedFile,'r') as F:
            for line in F:
                scaffoldID,start,stop,synthaseModelID,score,strand = line.strip().split('\t')
                assemblyID,chromID = scaffoldID.split('.')
                # print(assemblyID)
                if specificAssemblyID == assemblyID:
                    OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffoldID,start,stop,synthaseModelID,score,strand))
        
                        
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <bed file> <specific assembly ID> <molecule ID>\n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

bedFile = sys.argv[1]
specificAssemblyID = sys.argv[2]
moleculeID = sys.argv[3]

readBedFile(bedFile,specificAssemblyID,moleculeID)
