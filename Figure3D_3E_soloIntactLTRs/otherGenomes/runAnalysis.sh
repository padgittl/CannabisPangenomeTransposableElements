samtools faidx Ad77271a.a05.genome.fasta
cut -f 1,2 Ad77271a.a05.genome.fasta.fai > Ad77271a.chromSizes.txt

# identify solo LTRs for full assembly
python scripts/identifyScaffoldAndGlobalSoloLTRs.py Ad77271a.a05.genome.fasta.mod.EDTA.TEanno.gff3 Ad77271a.a05.genes.gff3 Ad77271a.chromSizes.txt 10 Ad77271a other 5000
## expected output files -->
## genomeID_Ty1_bothStrands_binWidth10_soloLTRFlankingWindow5000_candidateSoloLTRData.tsv
## genomeID_Ty1_bothStrands_binWidth10_soloLTRFlankingWindow5000_globalCandidateSoloLTRData.tsv
## genomeID_Ty3_bothStrands_binWidth10_soloLTRFlankingWindow5000_candidateSoloLTRData.tsv
## genomeID_Ty3_bothStrands_binWidth10_soloLTRFlankingWindow5000_globalCandidateSoloLTRData.tsv
## genomeID_unfilteredSoloLTRs.tsv

## *globalCandidateSoloLTRData.tsv is an intermediate file that gives the total count of solo-LTRs for the genome

# concatenate output files
cat *soloLTRFlankingWindow5000_candidateSoloLTRData.tsv > combinedFiltered_window5000bp.tsv
cat *unfilteredSoloLTRs.tsv > combinedUnfiltered.tsv

ls *TEanno*gff3 > gffFileList.txt

# calculate solo:intact ratio for scaffolds and globally
python scripts/calculateSoloIntactRatio.py combinedUnfiltered.tsv combinedFiltered_window5000bp.tsv gffFileList.txt 5000 other
## expected output files -->
## Ty3_globalSoloIntactRatio.tsv # full genome
## Ty1_globalSoloIntactRatio.tsv # full genome
## Ty1_soloIntactRatio.tsv # scaffold-level
## Ty3_soloIntactRatio.tsv # scaffold level
