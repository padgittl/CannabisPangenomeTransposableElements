# create conda environment
create_env.sh

# download gff files -- these include gene models and EDTA output
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/soloIntactLTRs/SODLb.primary_high_confidence.gff3.tar.gz
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/soloIntactLTRs/SODLb.unmasked.fasta.mod.EDTA.TEanno.gff3.tar.gz

tar -xvzf SODLb.primary_high_confidence.gff3.tar.gz
tar -xvzf SODLb.unmasked.fasta.mod.EDTA.TEanno.gff3.tar.gz

# run analysis to identify solo LTRs -- this script takes a while to run, so if testing in a hurry, can use chr7 only
# python scripts/identifySoloLTRs.py SODLb.unmasked.fasta.mod.EDTA.TEanno.gff3 SODLb.primary_high_confidence.gff3 10 SODLb cannabis 5000
python scripts/identifySoloLTRs.py data/SODLb.unmasked.fasta.mod.EDTA.TEanno.chr7.gff3 data/SODLb.primary_high_confidence.chr7.gff3 10 SODLb cannabis 5000
python scripts/identifySoloLTRs.py data/AH3Ma.unmasked.fasta.mod.EDTA.TEanno.chr1.gff3 data/AH3Ma.primary_high_confidence.chr1.gff3 10 AH3Ma cannabis 5000
python scripts/identifySoloLTRs.py data/AH3Mb.unmasked.fasta.mod.EDTA.TEanno.chr5.gff3 data/AH3Mb.primary_high_confidence.chr5.gff3 10 AH3Mb cannabis 5000

# concatenate output files
cat *soloLTRFlankingWindow5000_candidateSoloLTRData.tsv > combinedFiltered_window5000bp.tsv
cat *unfilteredSoloLTRs.tsv > combinedUnfiltered.tsv
# ls data/SODLb.unmasked.fasta.mod.EDTA.TEanno.chr7.gff3 > gffFileList.txt
ls data/*TEanno*gff3 > gffFileList.txt

# calculate solo:intact LTR ratio and generate boxplot per chromosome
python scripts/calculateSoloIntactRatio.py combinedUnfiltered.tsv combinedFiltered_window5000bp.tsv gffFileList.txt 5000 cannabis

# calculate solo:intact LTR ratio for full assembly
python scripts/calculateGlobalSoloIntactRatio.py combinedUnfiltered.tsv combinedFiltered_window5000bp.tsv gffFileList.txt 5000 cannabis
