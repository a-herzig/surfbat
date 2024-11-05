#!/bin/bash 
# SURFBAT 2024
##AF. Herzig




tmpdir=./
DIR=${3:-$tmpdir}
out=./results
OUT=${2:-$out}

x1=0.01
x2=0.7

X1=${4:-$x1}
X2=${5:-$x2}

label=${6}

invcf=$1

echo "CHROM POS ID REF ALT AC AN AF INFO beta varBeta chisq p" > ${OUT}.surfbat
echo "Extracting dosages"

zcat ${invcf} | grep -v "^#" | grep "SAP" | cut -f 1-5,8 --output-delimiter=" " | sed --regexp-extended 's/IMP;//g;s/AF=//g;s/INFO=//g;s/\;/ /g;s/AC=//g;s/AN=//g' | cut -d" " -f1-9 > ${tmpdir}SAP.${label}info.temp

k1=`zcat ${invcf} | grep -m1 "^#CHROM" | wc | awk '{print $2}'`
k1=$((5*($k1-9)))
seq 3 5 $k1 > ${tmpdir}cols.${label}.txt
seq 4 5 $k1 >> ${tmpdir}cols.${label}.txt

zcat ${invcf} | grep -v "^#" | grep "SAP" | cut -f 10- --output-delimiter=" " | tr "[ ]" "[:]" | cut -d":" -f$(echo `sort -n ${tmpdir}cols.${label}.txt` | tr "[ ]" "[,]") | tr "[,:]" "[  ]" > ${tmpdir}SAP.${label}.temp

rm ${tmpdir}cols.${label}.txt

echo "Surfbat test"

awk -v z1=${X1} -v z2=${X2} -v indexfile=${tmpdir}SAP.${label}info.temp 'BEGIN {j=1;while((getline < indexfile) >0) {maf[j]=$8;info[j]=$9;j++;}; j=1} {if(maf[j]>z1 && info[j]>z2) {{for(i=1;i<=NF/4;i++)(r1=r1+$(4*i-3)*(1-$(4*i-1))+$(4*i-2)*(1-$(4*i))) (r2=r2+$(4*i)*(1-$(4*i-2))+$(4*i-1)*(1-$(4*i-3))); if(r1>0 && r2>0){print log(r1/r2),(r2+r1)/(r1*r2),log(r1/r2)*log(r1/r2)/((r2+r1)/(r1*r2))} else {print "NA NA NA"};r1=0;r2=0}} else {print "NA NA NA"};j++}' SAP.${label}.temp > ${tmpdir}SAP.${label}stats

rm ${tmpdir}SAP.${label}.temp

Rscript surfbat_pvalue.R ${tmpdir}SAP.${label}stats ${tmpdir}SAP.${label}p
echo "Tidying"

paste -d' ' ${tmpdir}SAP.${label}info.temp ${tmpdir}SAP.${label}stats ${tmpdir}SAP.${label}p > ${tmpdir}SAP.${label}F

cat ${tmpdir}SAP.${label}F >> ${OUT}.surfbat

rm ${tmpdir}SAP.${label}*

echo "Fin"
