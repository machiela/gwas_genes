wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/technical/reference/genetic_map_b36.tar.gz
tar xvzf genetic_map_b36.tar.gz
gunzip ./genetic_map_b36/*.txt.gz
mv ./genetic_map_b36/*_b36.txt ./
rm genetic_map_b36.tar.gz
rm -rf genetic_map_b36
