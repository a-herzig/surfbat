#!/usr/bin/env Rscript
param = commandArgs(trailingOnly=TRUE)
options(warn=-1)
infile=param[1]
outfile=param[2]
write(pchisq(as.numeric(system(paste("awk '{print $3}' ",infile,sep=""),intern=T)),1,lower.tail=F),outfile,ncolumns=1)
