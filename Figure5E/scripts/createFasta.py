import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Abv,Chemotype,Type
# 79X,1,mj
def readChemotypeFile(chemotypeFile):
    populationDict = {}
    cleanPopIDs = {'mj':'mj', 'hc_hemp':'hc_hemp',
                   'feral':'feral', 'hemp':'hemp', 'asian_hemp':'Asian_hemp', 'F1':'F1'}
    with open(chemotypeFile,'r') as F:
        for line in F:
            if 'Chemotype' not in line and 'EXCLUDED' not in line and ',,' not in line:
                assemblyID,chemotypeID,populationID = line.strip().split(',')
                chemotypeID = 'Type' + chemotypeID
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


def readFastaFiles(dataList,synthaseData,specificChromID):
    seqDict = {}
    for fileID in dataList:
        for record in SeqIO.parse(fileID,"fasta"):
            if record.id in synthaseData:
                if record.id not in seqDict:
                    seqDict[record.id] = record
    return(seqDict)


def readBedToolsFile(bedToolsData,populationDict):
    synthaseData = {}
    lengthDict = {}
    teLengthThreshold = 100
    OUT_TE_lengths = open('te_lengths.txt','w')
    for fileID in bedToolsData:
        with open(fileID,'r') as F:
            for line in F:
                synthaseScaffoldID,synthaseStart,synthaseStop,synthaseID,column,synthaseStrand,teScaffoldID,teSource,teFeature,teStart,teStop,teScore,teStrand,teFrame,teAttribute,overlap = line.strip().split("\t")
                teStart = int(teStart)
                teStop = int(teStop)
                teLen = teStop - teStart + 1
                assemblyID,chrID = synthaseScaffoldID.split('.')
                if teLen >= teLengthThreshold:
                    if teFeature != 'repeat_region' and teFeature != 'long_terminal_repeat' and teFeature != 'target_site_duplication':
                        if 'LTR_retrotransposon' in teFeature and 'Gypsy' not in teFeature and 'Copia' not in teFeature:
                            teFeature = teFeature.replace('LTR_retrotransposon','unknown_LTR_retrotransposon')
                        if 'Gypsy' in teFeature:
                            teFeature = teFeature.replace('Gypsy','Ty3')
                        if 'Copia' in teFeature:
                            teFeature = teFeature.replace('Copia','Ty1')
                        if assemblyID in populationDict:
                            chemotypeID,popID = populationDict[assemblyID]
                        else:
                            print(assemblyID,"not in populationDict")
                            sys.exit()
                        chromID,moleculeID,countID = synthaseID.split('_')
                        attrCols = teAttribute.split(';')
                        #print(len(attrCols))
                        # MM3v1a ['ID=TE_struc_1928', 'Name=TE_00000016', 'Classification=DNA/DTC', 'Sequence_ontology=SO:0002285', 'Identity=1', 'Method=structural', 'TSD=ATT_ATT_100.0', 'TIR=CACTACAAGA_TCTTGTAGTG']
                        #print(assemblyID,attrCols)
                        # ['ID=LTRRT_9595', 'Parent=repeat_region_9595', 'Name=TE_00001807', 'Classification=LTR/Gypsy', 'Sequence_ontology=SO:0002265', 'ltr_identity=0.9973', 'Method=structural', 'motif=TGCA', 'tsd=TAGGA']
                        if len(attrCols) == 9:
                            getTEID = re.search('ID=(.+)',attrCols[0])
                            teID = getTEID.group(1)
                            getName = re.search('Name=(.+)',attrCols[2])
                            teName = getName.group(1)
                        else:
                            getTEID = re.search('ID=(.+)',attrCols[0])
                            teID = getTEID.group(1)
                            getName = re.search('Name=(.+)',attrCols[1])
                            teName = getName.group(1)
                            #print(teID,teName)
                        # newTEID = teFeature + '_' + teID + '_' + teName + '_' + synthaseID + '_' + chemotypeID + '_' + popID
                        # newTEID = teID + '_' + teFeature + '_' + teName + '_' + synthaseID + '_' + chemotypeID + '_' + popID
                        teFeature = teFeature.replace('_retrotransposon','')
                        teID = teID.replace('TE_homo_','')
                        newTEID = teID + '_' + synthaseID + '_' + teFeature + '_' + chemotypeID + '_' + popID
                        if synthaseScaffoldID not in synthaseData:
                            synthaseData[synthaseScaffoldID] = {}
                        if moleculeID not in synthaseData[synthaseScaffoldID]:
                            synthaseData[synthaseScaffoldID][moleculeID] = {}
                        if teFeature not in synthaseData[synthaseScaffoldID][moleculeID]:
                            synthaseData[synthaseScaffoldID][moleculeID][teFeature] = {}
                        if teID not in synthaseData[synthaseScaffoldID][moleculeID][teFeature]:
                            synthaseData[synthaseScaffoldID][moleculeID][teFeature][teID] = []
                        synthaseData[synthaseScaffoldID][moleculeID][teFeature][teID].append((synthaseID,synthaseStrand,newTEID,teStart,teStop))
                        if moleculeID not in lengthDict:
                            lengthDict[moleculeID] = {}
                        if teFeature not in lengthDict[moleculeID]:
                            lengthDict[moleculeID][teFeature] = []
                        lengthDict[moleculeID][teFeature].append(teLen)
                        #print(teLen)
                        OUT_TE_lengths.write("%s\t%s\t%s\n" % (moleculeID,teFeature,teLen))
    return(synthaseData,lengthDict)


def createFasta(seqDict,synthaseData):
    catCommands = {}
    commandDict = {}
    # synthaseData[synthaseScaffoldID][moleculeID][teFeature][teID].append((synthaseID,synthaseStrand,newTEID,teStart,teStop))
    for synthaseScaffoldID in synthaseData:
        if synthaseScaffoldID in seqDict:
            record = seqDict[synthaseScaffoldID]
            for	moleculeID in synthaseData[synthaseScaffoldID]:
                if moleculeID not in commandDict:
                    commandDict[moleculeID] = {}
                for teFeature in synthaseData[synthaseScaffoldID][moleculeID]:
                    OUT = synthaseScaffoldID + '_' + moleculeID + '_' + teFeature + '.fasta'
                    cat_command = "cat " + OUT + " >> " + moleculeID + "_" + teFeature + ".fasta"
                    if moleculeID not in catCommands:
                        catCommands[moleculeID] = {}
                    if teFeature not in catCommands[moleculeID]:
                        catCommands[moleculeID][teFeature] = []
                    catCommands[moleculeID][teFeature].append(cat_command)
                    mafft_command = "mafft --auto " + moleculeID + "_" + teFeature + ".fasta > " + moleculeID + "_" + teFeature + "_aln.fasta"
                    # clustalw2 -infile=CBCAS_unknown_LTR_retrotransposon.fasta -type=DNA -outfile=CBCAS_unknown_LTR_retrotransposon.clustal -output=NEXUS
                    # https://bioinformaticsworkbook.org/phylogenetics/FastTree.html#gsc.tab=0
                    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894340/
                    # "Alternative settings for FastTree include JC (the Jukes-Cantor model) [36] instead of GTR, but this simplified model is not recommended except under very unusual circumstances where the data seem to fit the Jukes-Cantor model best (unlikely for most data). Note that the GAMMA setting is usually used in phylogenetic analyses"
                    # FastTree -nt -gtr -gamma input_aln.fasta > input.tre
                    fasttree_command = "FastTree -nt -gtr -gamma " + moleculeID + "_" + teFeature + "_aln.fasta > " + moleculeID + "_" + teFeature + "_aln.tree" 
                    if teFeature not in	commandDict[moleculeID]:
                        commandDict[moleculeID][teFeature] = (mafft_command,fasttree_command)
                    recordDict = {}
                    recordList = []
                    for teID in synthaseData[synthaseScaffoldID][moleculeID][teFeature]:
                        if len(synthaseData[synthaseScaffoldID][moleculeID][teFeature][teID]) > 1:
                            synthaseID,synthaseStrand,seqLabel,teStart,teStop = synthaseData[synthaseScaffoldID][moleculeID][teFeature][teID][0]
                            subSeq = record.seq[teStart-1:teStop]
                            if synthaseStrand == '+':
                                sequence = str(subSeq)
                            else:
                                sequence = subSeq.reverse_complement()
                                sequence = str(sequence)
                            newRecord = SeqRecord(Seq(sequence), id=seqLabel, name='', description='')
                            if seqLabel not in recordDict:
                                recordDict[seqLabel] = newRecord
                                recordList.append(newRecord)
                        else:
                            for synthaseID,synthaseStrand,seqLabel,teStart,teStop in synthaseData[synthaseScaffoldID][moleculeID][teFeature][teID]:
                                # record.seq[1-1:len(record.seq)]
                                subSeq = record.seq[teStart-1:teStop]
                                if synthaseStrand == '+':
                                    sequence = str(subSeq)
                                else:
                                    sequence = subSeq.reverse_complement()
                                    sequence = str(sequence)
                                    # record = SeqRecord(Seq("ATGCCCGGGAAATAG"))
                                    # SeqRecord(seq=Seq('ATGCCCGGGAAATAG'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])
                                newRecord = SeqRecord(Seq(sequence), id=seqLabel, name='', description='')
                                if seqLabel not in recordDict:
                                    recordDict[seqLabel] = newRecord
                                    recordList.append(newRecord)
                    SeqIO.write(recordList,OUT,"fasta")
                            
                    
        else:
            print(synthaseScaffoldID,'not in seqDict')
            sys.exit()
    # commandDict[moleculeID][teFeature] = (cat_command,mafft_command,fasttree_command)
    for moleculeID in commandDict:
        for teFeature in commandDict[moleculeID]:
            OUT_COMMANDS = open(moleculeID + "_" + teFeature + "_run_commands.sh",'w')
            mafft_command,fasttree_command = commandDict[moleculeID][teFeature]
            OUT_COMMANDS.write("%s\n" % (mafft_command))
            OUT_COMMANDS.write("%s\n" % (fasttree_command))
    # catCommands[moleculeID][teFeature].append(cat_command)
    for moleculeID in catCommands:
        for teFeature in catCommands[moleculeID]:
            CAT_OUT = open(moleculeID + "_" + teFeature + '_cat_command.sh','w')
            for cat_command in catCommands[moleculeID][teFeature]:
                CAT_OUT.write("%s\n" % (cat_command))


            
'''
record.seq = record.seq[:60]
print record.id,record.seq
if record.id not in recordDict:
recordDict[record.id] = 1
recordList.append(record)

for scaffoldID in seqDict:
seqRecord = seqDict[scaffoldID]
'''        

# lengthDict[moleculeID][teFeature].append(teLen)
def createHist(lengthDict):
    binWidth = 500
    for moleculeID in lengthDict:
        for teFeature in lengthDict[moleculeID]:
            #if len(lengthDict[moleculeID][teFeature]) > 2:
            #dataList = []
            #print(len(lengthDict[moleculeID][teFeature]),lengthDict[moleculeID][teFeature])
            #for length in lengthDict[moleculeID][teFeature]:
            #print(moleculeID,teFeature,length)
                #dataList.append(length)
                #print(dataList)
                '''
                fig = plt.figure(figsize=(5,4))
                fullList = dataList
                #print(dataList)
                bins = np.arange(min(fullList),max(fullList),binWidth)
                counts1, bins1, bars1 = plt.hist(fullList, bins = bins, histtype = 'stepfilled', color ='#0471A6', label = moleculeID + ' ' + teFeature, alpha=0.75, linewidth=2)
                plt.xlabel('Transposable element length (bp)',size=16)
                plt.ylabel('Count',size=16)
                plt.legend(frameon=False,fontsize=14)
                plt.tight_layout()
                plt.savefig(moleculeID + '_' + teFeature + '_lengths.png' , dpi=600)
                plt.close()
                '''
                
########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <fasta file list> <bedtools file list> <chemotype/population ID file> <specific chromosome ID>\n"
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

fastaFileList = sys.argv[1]
bedToolsFileList = sys.argv[2]
chemotypeFile = sys.argv[3]
specificChromID = sys.argv[4]

populationDict = readChemotypeFile(chemotypeFile)

dataList = readFileList(fastaFileList)

bedToolsData = readFileList(bedToolsFileList)
synthaseData,lengthDict = readBedToolsFile(bedToolsData,populationDict)

seqDict = readFastaFiles(dataList,synthaseData,specificChromID)
createFasta(seqDict,synthaseData)
#createHist(lengthDict)
