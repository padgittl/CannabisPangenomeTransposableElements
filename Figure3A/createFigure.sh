gunzip data/*gz
ls data/*sum > sumFileList.txt

# generate csv file for plotting
python scripts/createCSVForJoyplot.py sumFileList.txt ../utilityFiles/pangenome_chemotype_population_IDs.csv

# plot TE percentages as a joyplot
python scripts/joyplot.py dnaTE_percentMasked.csv popID_teID
