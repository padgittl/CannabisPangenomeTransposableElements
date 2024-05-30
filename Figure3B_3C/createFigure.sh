# create conda environment
create_env.sh

# create file list
ls data/*gff3 > gffFileList.txt

# create figure of TE time divergence
python scripts/timeDistribution.py gffFileList.txt
