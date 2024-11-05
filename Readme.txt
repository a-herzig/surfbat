###Readme For Surfbat
##Contact anthony.herzig@inserm.fr for more info

In this directory you can find:

1) A markdown (Surfbat.example.md) explaining the msprime simulations from the original surfbat article, going from data creating until running the association testdata/5popSim_4B11_Cases_test
2) Test data for quickly testing surfbat in /testdata
You can use the following command once Impute5 has been installed:

impute5_static	\
	--h testdata/5popSim_4B11_Ref.bcf \
	--g testdata/5popSim_4B11_Cases_test.bcf \
	--r 1:1-102000000 \
	--buffer-region 1:1-102000000 \
	--surfbat --o testdata/5popSim_4B11_Cases_test_out.vcf
	
3) Summary statistics from Figure 2 in the original article
4) Bash and R scripts for calculating the surfbat test statistics by hand, details within the example markdown