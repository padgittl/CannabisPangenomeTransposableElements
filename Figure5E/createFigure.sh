# create conda environment
create_env.sh

# create genome-centric bed files for the cannabinoid synthase annotations
less utilityFiles/genomeIDs.txt | awk '{print "python scripts/createSeparateBedFiles.py utilityFiles/CBCAS_fullhits.bed "$1" CBCAS"}' > getSeparateCBCAS.sh
less utilityFiles/genomeIDs.txt | awk '{print "python scripts/createSeparateBedFiles.py utilityFiles/CBDAS_fullhits.bed "$1" CBDAS"}' > getSeparateCBDAS.sh
less utilityFiles/genomeIDs.txt | awk '{print "python scripts/createSeparateBedFiles.py utilityFiles/THCAS_fullhits.bed "$1" THCAS"}' > getSeparateTHCAS.sh

# download data
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb.softmasked.chr7.fasta
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb.unmasked.fasta.mod.EDTA.TEanno.chr7.gff3
wget https://salk-tm-pub.s3.us-west-2.amazonaws.com/cannabis_pangenome/te_landscape_analysis/SODLb.unmasked.fasta.mod.EDTA.TEanno.gff3.tar.gz
tar -xzvf SODLb.unmasked.fasta.mod.EDTA.TEanno.gff3.tar.gz 

# download genome and other data from figshare
https://figshare.com/articles/dataset/SODLb_tar_gz/25877290

# create file of chromosome lengths
for f in *.softmasked.fasta; do echo samtools faidx $f; done | grep "a\." > createChromIndex_a.sh
for f in *.softmasked.fasta; do echo samtools faidx $f; done | grep "b\." > createChromIndex_b.sh
for f in *.softmasked.fasta.fai; do echo cut -f 1,2 $f '>' `basename $f .softmasked.fasta.fai`.chromSizes.txt; done > getChromSizes.sh

# retrieve genomic windows that are flanking cannabinoid synthases with bedtools flank
for f in *_CBCAS_fullhits.bed; do echo bedtools flank -i $f -g `basename $f _CBCAS_fullhits.bed`.chromSizes.txt -l 2000 -r 2000 '>' `basename $f _CBCAS_fullhits.bed`_CBCAS_flanking_2000.bed; done > flank_CBCAS.sh
for f in *_CBDAS_fullhits.bed; do echo bedtools flank -i $f -g `basename $f _CBDAS_fullhits.bed`.chromSizes.txt -l 2000 -r 2000 '>' `basename $f _CBDAS_fullhits.bed`_CBDAS_flanking_2000.bed; done > flank_CBDAS.sh
for f in *_THCAS_fullhits.bed; do echo bedtools flank -i $f -g `basename $f _THCAS_fullhits.bed`.chromSizes.txt -l 2000 -r 2000 '>' `basename $f _THCAS_fullhits.bed`_THCAS_flanking_2000.bed; done > flank_THCAS.sh

# prepare to get TEs that overlap with cannabinoid synthases with bedtools intersect
for f in *_CBCAS_flanking_2000.bed; do echo bedtools intersect -a $f -b `basename $f _CBCAS_flanking_2000.bed`.unmasked.fasta.mod.EDTA.TEanno.gff3 -wo '>' `basename $f _CBCAS_flanking_2000.bed`_CBCAS_intersect_2000.bed; done > intersect_CBCAS.sh
for f in *_CBDAS_flanking_2000.bed; do echo bedtools intersect -a $f -b `basename $f _CBDAS_flanking_2000.bed`.unmasked.fasta.mod.EDTA.TEanno.gff3 -wo '>' `basename $f _CBDAS_flanking_2000.bed`_CBDAS_intersect_2000.bed; done > intersect_CBDAS.sh
for f in *_THCAS_flanking_2000.bed; do echo bedtools intersect -a $f -b `basename $f _THCAS_flanking_2000.bed`.unmasked.fasta.mod.EDTA.TEanno.gff3 -wo '>' `basename $f _THCAS_flanking_2000.bed`_THCAS_intersect_2000.bed; done > intersect_THCAS.sh

# run bedtools intersect
nohup ./intersect_CBCAS.sh > intersect_CBCAS.log 2>&1 &
nohup ./intersect_CBDAS.sh > intersect_CBDAS.log 2>&1 &
nohup ./intersect_THCAS.sh > intersect_THCAS.log 2>&1 &

# create list of genome fasta files
for f in *_CBCAS_flanking_2000.bed; do echo `basename $f _CBCAS_flanking_2000.bed`.softmasked.fasta; done > fastaFileList.txt
for f in *_CBDAS_flanking_2000.bed; do echo `basename $f _CBDAS_flanking_2000.bed`.softmasked.fasta; done >> fastaFileList.txt
for f in *_THCAS_flanking_2000.bed; do echo `basename $f _THCAS_flanking_2000.bed`.softmasked.fasta; done >> fastaFileList.txt
less fastaFileList.txt | sort | uniq > sorted_fastaFileList.txt

# create list of bedtools intersect output files
ls *_intersect_2000.bed > intersectFileList.txt

# retrieve genomic sequence with custom script
nohup python scripts/createFasta.py sorted_fastaFileList.txt intersectFileList.txt ../utilityFiles/pangenome_chemotype_population_IDs.csv chr7 > createFasta.log 2>&1 &
chmod +x *sh
for f in *_cat_command.sh; do ./$f; done

# create alignment with mafft and tree with FastTree
for f in *_run_commands.sh; do echo nohup ./$f ">" `basename $f .sh`.log "2>&1 &" ; done > run.sh
