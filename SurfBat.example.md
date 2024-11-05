```bash
##Here we use stdpopsim to simulate a tree sequence from the American Adixture scenario
stdpopsim HomSap -c 15 -o afr-america-chr15.trees -s 1213 -g HapMapII_GRCh38 -d AmericanAdmixture_4B11 AFR:15000 EUR:15000 ASIA:15000 ADMIX:2000
```

```python
import sys
import msprime 
import numpy as np
import stdpopsim
from matplotlib import pyplot as plt
import tskit
maf=0.05
##Read in the simulated tree sequence
vcf_name = "afr.america.chr15.msp"
ts = tskit.load("afr-america-chr15.trees")
model = msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES)

site_to_mask_index = [False] * ts.get_num_sites()
##Here we both remove rare variants and also remove every 5th variant in order to have simulated data that is not too big
for variant in ts.variants():
	freq = variant.frequencies()[next(iter(variant.frequencies()))]
	if (freq < maf) or (freq > 1-maf) or (variant.index % 5 > 0):
		site_to_mask_index[variant.index]=True

with open("".join([vcf_name,".vcf"]), "w") as vcf_file:
	ts.write_vcf(vcf_file, site_mask = site_to_mask_index)

quit()

```

```bash
##Probably this can be done in python, but we noticed that the ts.write_vcf command also puts the data on chromosme 1, so we change that by hand to chromsome 15
##We also remove multi-alleleic sites by hand
chr=15
sed -i "s/ID=1,/ID=${chr},/g" afr.america.chr${chr}.msp.vcf
awk -v chr=${chr} 'OFS="\t" {gsub(/^1/,chr); print}' afr.america.chr${chr}.msp.vcf | awk '/^#/ || !/,/' | bgzip > afr.america.chr${chr}.msp.vcf.gz
#rm afr.america.chr${chr}.msp.vcf
#rm afr-america-chr${chr}.trees
tabix afr.america.chr${chr}.msp.vcf.gz


```

```R
##Now we bring the data into R and treat it with gaston
##First step, restrict to very common variants, prune the data, randomly select 5,000 Snps to serve as a genotyping array 
library(gaston)
chr<-15
gwas<-read.vcf("afr.america.chr15.msp.vcf.gz")
chip<-gwas[,which(gwas@snps$maf>0.1)]
chip<-LD.thin(chip,0.8)
chip<-chip[,sort(sample(1:ncol(chip),5000,F))]


### What proceeds are sets of commands to take the total simulated data and make new vcf files, correctly formatted for IMPUTE5 that represent a Target dataset with only the array ###genotyping positions for the 2,000 admixed individuals, and the reference panel with all SNPs for the other 45,000 individuals from AFR, ASIA, and EUR.

###requirements, bcftools, shapeit5 and xcftools
write.table(cbind(15,chip@snps$pos),"chip.txt",col.names=F,row.names=F,sep="\t",quote=F)
gwas@ped$population<-c(rep("AFR",15000),rep("EUR",15000),rep("ASIA",15000),rep("ADMIX",2000))

A1<-gwas@ped$id[which(gwas@ped$population%in%c("AFR","EUR","ASIA"))]
A2<-gwas@ped$id[which(!gwas@ped$population%in%c("AFR","EUR","ASIA"))]

write.table(A1,"ref.txt",col.names=F,row.names=F,sep="\t",quote=F)
write.table(A2,"case.txt",col.names=F,row.names=F,sep="\t",quote=F)


system(paste("bcftools view -S ref.txt afr.america.chr15.msp.vcf.gz | bgzip > afr.america.chr15.msp_Ref.vcf.gz",sep=""))
system(paste("tabix -f afr.america.chr15.msp_Ref.vcf.gz",sep=""))


system(paste("bcftools view -S case.txt -T chip.txt afr.america.chr15.msp.vcf.gz | bgzip > afr.america.chr15.msp_Cases.vcf.gz",sep=""))
system(paste("tabix -f afr.america.chr15.msp_Cases.vcf.gz",sep=""))

system(paste("bcftools convert -O b -o afr.america.chr15.msp_Cases.bcf afr.america.chr15.msp_Cases.vcf.gz",sep=""))
system(paste("tabix -f afr.america.chr15.msp_Cases.bcf",sep=""))


system(paste("bcftools +fill-tags afr.america.chr15.msp_Cases.bcf -Ob -o temp.bcf -- -t AN,AC",sep=""))
system(paste("mv temp.bcf afr.america.chr15.msp_Cases.bcf",sep=""))
system(paste("tabix -f afr.america.chr15.msp_Cases.bcf",sep=""))

###We need the xcftools functionality of shapeit5
system(paste("/shapeit5/xcftools/bin/xcftools view --input afr.america.chr15.msp_Cases.bcf -O sh -o afr.america.chr15.msp_Cases_xcf.bcf --r ",15," --m 0",sep=""))

system(paste("bcftools convert -O b -o afr.america.chr15.msp_Ref.bcf afr.america.chr15.msp_Ref.vcf.gz",sep=""))
system(paste("tabix -f afr.america.chr15.msp_Ref.bcf",sep=""))

system(paste("bcftools +fill-tags afr.america.chr15.msp_Ref.bcf -Ob -o temp.bcf -- -t AN,AC",sep=""))
system(paste("mv temp.bcf afr.america.chr15.msp_Ref.bcf",sep=""))
system(paste("tabix -f afr.america.chr15.msp_Ref.bcf",sep=""))

system(paste("/shapeit5/xcftools/bin/xcftools view --input afr.america.chr15.msp_Ref.bcf -O sh -o afr.america.chr15.msp_Ref_xcf.bcf --r ",15," --m 0",sep=""))

ref1<-paste0("afr.america.chr15.msp_Ref_xcf")
target<-paste("afr.america.chr15.msp_Cases",sep="")

###Finally we can start to run IMPUTE5, firstly use the chunking algorithm
system(paste("imp5Chunker_v1.2.0_static --h ",ref1,".bcf --g ",target,"_xcf.bcf --r ",chr," --o coordinates.I51.2_",chr,".txt",sep=""))
####
###
##
####
###
##
corrdinates<-read.table(paste("coordinates.I51.2_",chr,".txt",sep=""),header=F,as.is=T)

###Ru IMPUTE5 over each chunk
for (j in 1:nrow(corrdinates)){

system(paste("impute5_static --h ",ref1,".bcf --g ",target,"_xcf.bcf --r ",corrdinates[j,4]," --m /PUBLIC_DATA/ReferencePanels/1kG/1000GP_Phase3_haplotypes/genetic_map_chr",chr,"_combined_b37.txt --surfbat --buffer-region ",corrdinates[j,3]," --o ",target,".I51.2_chunk",corrdinates[j,1],".vcf.gz --l ",target,".I51.2_chunk",corrdinates[j,1],".log",sep=""))

if(j==1){
write(paste(target,".I51.2_chunk",corrdinates[j,1],".vcf.gz",sep=""),paste("ligate.I51.2_chr",chr,".txt",sep=""))
} else {
write(paste(target,".I51.2_chunk",corrdinates[j,1],".vcf.gz",sep=""),paste("ligate.I51.2_chr",chr,".txt",sep=""),append=T)

}
}

### a bit or merging and tidying
system(paste("bcftools concat -n -fligate.I51.2_chr",chr,".txt -Oz -o ",target,".I51.2_chunkAll.vcf.gz",sep=""))

system(paste("tabix -f ",target,".I51.2_chunkAll.vcf.gz",sep=""))

for (j in 1:nrow(corrdinates)){
system(paste("rm ",target,".I51.2_chunk",corrdinates[j,1],".vcf.*",sep=""))
system(paste("rm ",target,".I51.2_chunk",corrdinates[j,1],".log",sep=""))
} 

system(paste("rm ligate.I51.2_chr",chr,".txt",sep=""))
system(paste("rm coordinates.I51.2_",chr,".txt",sep=""))


###Impute5 can directly compute the p-values, but we also supply a bash script to do so, which can help explin how the calculation is done for those interested
system(paste("./surfbat.i512.sh ",target,".I51.2_chunkAll.vcf.gz ",target,".I51.2_chunkAll /tmp/ 0.01 0.3 chr15.msp",sep=""))

### Read in the data, caclulate lamdba for example
sf<-read.table(paste(target,".I51.2_chunkAll.surfbat",sep=""),header=T,as.is=F)
sf<-sf[-which(is.na(sf[,13])),]
chisq1 <- sf[,12]
lambda1 = median(chisq1)/qchisq(0.5,1)

###fin

```