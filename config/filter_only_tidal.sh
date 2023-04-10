# retains only strains analyzed in TIDAL

grep -f ../../tetf_tidal/results/overlaps/tidal_strains.tsv sample_table.csv > sample_table.tidal.csv
grep -f ../../tetf_tidal/results/overlaps/tidal_strains.tsv subsample_table.csv > subsample_table.tidal.csv    