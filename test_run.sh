# Download test data for AH3Mb chromosome 1
mkdir data
aws s3 cp --recursive s3://salk-tm-dev/cannabis/lillian/CannabisPangenomeTransposableElements/AH3Mb/testData/ data/.

python scripts/plotTEsAndMethylation.py data/AH3Mb_gene_densities.txt data/AH3Mb_TEs_densities.txt data/AH3Mb_TEs_times.txt data/AH3Mb_LTRs_densities.txt data/AH3Mb_LTRs_times.txt data/AH3Mb_methylation_densities.txt data/AH3Mb_gc_info.txt data/AH3Mb_dups_variant_density.txt data/AH3Mb_inv_variant_density.txt data/AH3Mb_transloc_variant_density.txt data/AH3Mb_chromosomeLengths.txt AH3Mb 1000000