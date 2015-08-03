# get_recomb_genes.py
#
# Function:
# Takes a file in format (chromosome \t location \t name) and retrieves surrounding genes using a recombination map
# Calculates centiMorgans from LD block
# Uses hg18 coordinates
#
# Variables:
# cancer - name of cancer; for finding the file and naming the output
# threshold - minimum recombination rate (cM/Mb) used to define boundaries of LD blocks
# plusminus - number of base pairs extended from the LD boundaries to define the region of interest 
import sys

cancer=str(sys.argv[1]).strip()
threshold=20 
plusminus=500000

map_dir="/recomb/" # 1000 Genomes CEU recombination map directory (hg18)
gene_dir="/genes/" # directory for RefSeq genes divided by chromosome (hg18)
gwas_snps=open(cancer+"_snps.txt").readlines() # file containing significant GWAS SNPs (hg18)
out=open(cancer+"_ld.txt","w") # output file of LD regions (hg18)

header=["cancer","rsnum","chr","snp_bp","snp_cM","ld_block_start","ld_block_end","ld_window_start","ld_window_end","gene","transcript_id","tx_start","tx_end","ld_block","gene_bp","gene_cM","gene_bp_ld"]
print >> out,"\t".join(header)

for snp in gwas_snps:
	rsnum,chr,loc=snp.strip().split()
	loc=int(loc)
	
	# forms the recombination map only including the base pairs that exceed the threshold value (cM/Mb)
	recomb_map=open(map_dir+'genetic_map_'+chr+'_b36.txt').readlines() # opens the recombination map file
	recomb_dict={} # initializes dictionary containing base pair number and recombination rate
	recomb_peaks=[] # initializes list for bp with cM/Mb exceeding threshold
	
	for j in range(1,len(recomb_map)):
		bp,rate,cM=recomb_map[j].strip().split()
		recomb_dict[bp]=float(cM)
		if float(rate)>=threshold:
			recomb_peaks.append(int(bp))
		elif j==1: # caps beginning of list with first entry of recombination map
			recomb_peaks.append(int(bp))
		elif j==len(recomb_map)-1: # caps the end of the list with the last entry in the recombination map
			recomb_peaks.append(int(bp))
	
	# Get base pair coordinates from recombination dictionary
	recomb_bp=sorted([int(i) for i in recomb_dict.keys()]) # list of all bp in ascending order
	
	# Find the LD block and window around the GWAS SNP
	k=0
	while k<len(recomb_peaks)-1 and recomb_peaks[k]<loc:
		k+=1
		
	
	ld_block=(recomb_peaks[k-1],recomb_peaks[k]) # sets the LD block
	ld_window=[ld_block[0]-plusminus,ld_block[1]+plusminus] 
	
	# Adjust bounds for LD window if they exceed the chromosome boundries
	if ld_window[0]<recomb_bp[0]: 
		ld_window[0]=recomb_bp[0]
	if ld_window[1]>recomb_bp[-1]:
		ld_window[1]=recomb_bp[-1]
	
	# Calculate the distance in centimorgans (cM) the SNP is from the LD block boundaries
	k=0
	while k<len(recomb_bp)-1 and recomb_bp[k]<loc: 
		k+=1
	
	closest_bp=min(loc-recomb_bp[k-1],recomb_bp[k]-loc)
	if loc-recomb_bp[k-1]==closest_bp:
		snp_cM=recomb_dict[str(recomb_bp[k-1])]
	else:
		snp_cM=recomb_dict[str(recomb_bp[k])]
	
	ld_block_cM=[recomb_dict[str(ld_block[1])],recomb_dict[str(ld_block[1])]]
	ld_block_cM_SNP=[ld_block_cM[0]-snp_cM,ld_block_cM[1]-snp_cM]
	
	
	# Find genes located within the LD window around the SNP
	genes=open(gene_dir+"refGene_"+chr+".txt").readlines() # opens refGene file of genes for chromosome
	for i in range(0,len(genes)): # loop goes through all genes in the file
		gene=genes[i].strip().split()
		txStart=int(gene[4])
		txEnd=int(gene[5])
		
		# Proceeds if gene spans the LD window
		if (ld_window[0]<txStart<ld_window[1]) or (ld_window[0]<txEnd<ld_window[1]):
			name=gene[1]
			name2=gene[12]
			
			# Calculate distance and centimorgans away from the GWAS SNP
			if (ld_block[0]<=txStart<=ld_block[1]) or (ld_block[0]<=txEnd<=ld_block[1]) or (txStart<=ld_block[0]<=ld_block[1]<=txEnd):
				in_ld_block="1"
				gene_bp_ld=0
			else:
				in_ld_block="0"
				if txStart>ld_block[1]:
					gene_bp_ld=txStart-ld_block[1]
				elif ld_block[0]>txEnd:
					gene_bp_ld=(ld_block[0]-txEnd)*(-1)
			
			if txStart<=loc<=txEnd:
				gene_bp=0
				gene_cM=0
			
			elif txEnd<=loc:
				gene_bp=txEnd-loc
				k=0
				while k<len(recomb_bp)-1 and recomb_bp[k]<txEnd: 
					k+=1
				closest_bp=min(txEnd-recomb_bp[k-1],recomb_bp[k]-txEnd)
				if txEnd-recomb_bp[k-1]==closest_bp:
					gene_cM=recomb_dict[str(recomb_bp[k-1])]-snp_cM
				else:
					gene_cM=recomb_dict[str(recomb_bp[k])]-snp_cM
			
			elif txStart>=loc:
				gene_bp=txStart-loc
				k=0
				while k<len(recomb_bp)-1 and recomb_bp[k]<txStart: 
					k+=1
				closest_bp=min(txStart-recomb_bp[k-1],recomb_bp[k]-txStart)
				if txStart-recomb_bp[k-1]==closest_bp:
					gene_cM=recomb_dict[str(recomb_bp[k-1])]-snp_cM
				else:
					gene_cM=recomb_dict[str(recomb_bp[k])]-snp_cM
			
			print_vars=[cancer,rsnum,chr,loc,snp_cM,ld_block[0],ld_block[1],ld_window[0],ld_window[1],name2,name,txStart,txEnd,in_ld_block,gene_bp,gene_cM,gene_bp_ld]
			print_vars_str=[str(i) for i in print_vars]
			print >> out,"\t".join(print_vars_str)


out.close()
