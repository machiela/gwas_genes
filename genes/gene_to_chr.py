# gene_to_chr.py
#
# This program divided the hg18 RefSeq file into smaller lists by chromosome.
# These smaller files are accessed by get_recomb_genes.py to find the genes within the target region.

import subprocess
subprocess.call("wget ftp://hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg18/database/refGene.txt.gz", shell=True)
subprocess.call("gunzip refGene.txt.gz", shell=True)

genes=open("refGene.txt").readlines()

for i in range(0,len(genes)):
	gene=genes[i].strip().split()
	chr=gene[2]
	chr_out=open("refGene_"+chr+".txt","a")
	print >> chr_out,genes[i].strip()
	chr_out.close()

