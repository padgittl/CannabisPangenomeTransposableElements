# create conda environment
create_env.sh

# download files
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb.loreme.chr7.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb.primary_high_confidence_chr7.gff3
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb.softmasked.chr7.fasta
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb.unmasked.fasta.mod.EDTA.TEanno.chr7.gff3
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb_dups_chr7.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb_invs_chr7.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb_trans_chr7.bed

# methylation data is available for download for the following genomes (each file is ~1 GB):
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/AH3Ma.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/AH3Mb.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/FCS1a.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/FCS1b.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/H3S1a.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/H3S1b.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/KCDv1a.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/KCDv1b.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/NLv1a.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/NLv1b.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/SAN2a.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/SAN2b.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/SN1v3a.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/SN1v3b.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/SODLa.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/SODLb.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/WHWa.loreme.bed
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/csat_methylation/WHWb.loreme.bed


# process data for figure -- this script takes a few minutes to run
python scripts/teAndMethylationDensities.py SODLb.primary_high_confidence_chr7.gff3 SODLb.unmasked.fasta.mod.EDTA.TEanno.chr7.gff3 SODLb.loreme.chr7.bed SODLb.softmasked.chr7.fasta data/other_files/csat_orientations.tsv SODLb_dups_chr7.bed SODLb_invs_chr7.bed SODLb_trans_chr7.bed SODLb chr7 1000000 &

# create figure -- this script takes a few minutes to run
python scripts/plotTELandscape.py data/postprocessed/SODLb_gene_densities_window1000000_chr7.txt data/postprocessed/SODLb_TEs_densities_window1000000_chr7.txt data/postprocessed/SODLb_TEs_times_window1000000_chr7.txt data/postprocessed/SODLb_LTRs_densities_window1000000_chr7.txt data/postprocessed/SODLb_LTRs_times_window1000000_chr7.txt data/postprocessed/SODLb_methylation_densities_window1000000_chr7.txt data/postprocessed/SODLb_gc_info_window1000000_chr7.txt data/postprocessed/SODLb_dups_variant_densities_window1000000_chr7.txt data/postprocessed/SODLb_inv_variant_densities_window1000000_chr7.txt data/postprocessed/SODLb_transloc_variant_densities_window1000000_chr7.txt data/postprocessed/SODLb_chromosomeLengths_window1000000_chr7.txt data/postprocessed/SODLb_Ty1_soloIntactRatio_window1000000_chr7.tsv data/postprocessed/SODLb_Ty3_soloIntactRatio_window1000000_chr7.tsv data/other_files/allSynthaseHits.bed data/other_files/pangenome.ALT4.labeled.bed data/other_files/pangenome.BKR.labeled.bed data/other_files/csat_orientations.tsv SODLb 1000000 0 70719495 PIF_Harbinger_TIR_transposon Ty1_LTR_retrotransposon
