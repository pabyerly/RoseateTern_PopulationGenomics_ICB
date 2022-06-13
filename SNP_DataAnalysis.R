library("adegenet")
library("hierfstat")
library("pegas")
library(vcfR)
library(reshape2)
library(ggplot2)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(dplyr)
library(diveRsity)
library(nlme)
library("PMCMRplus")
library(rstatix)

#roseate tern palette
myCol=c("deepskyblue3","grey44", "orangered","hotpink1", "black")
#################################################################################################

#adegenet: https://grunwaldlab.github.io/Population_Genetics_in_R/reading_vcf.html
#convert vcf to genlight via https://knausb.github.io/vcfR_documentation/export_genlight_snpclone.html
#input vcf file

vcf_mod<- read.vcfR("ROST_SNPS_ModerateFilterSet.vcf")

vcf_strict<- read.vcfR("ROST_SNPS_StrictFilterSet.vcf")

#to genind format
rost_vcf_mod<-vcfR2genind(vcf_mod)

rost_vcf_strict<-vcfR2genind(vcf_strict)

#read in population file 
snp_pop=read.csv("snp_pop.csv")

#assign population, ie put it in the population slot 
pop(rost_vcf_mod)=snp_pop05$Pop
pop(rost_vcf_strict)=snp_pop05$Pop

###############################################################################
#HWE test
#chi squared test of HWE, 2 p-values: 1 analytical, 1 from permutations
(snps_hwe<-hw.test(rost_vcf_mod, B=1000))
write.csv(x=snps_hwe, file = "mod_HWE")

(snps_hwe<-hw.test(rost_vcf_strict, B=1000))
write.csv(x=snps_hwe, file = "strict_HWE")

#check for each population
(snps_hwe_pop<-seppop(rost_vcf_mod) %>% lapply (hw.test, B = 1000))

(snps_hwe_pop<-seppop(rost_vcf_strict) %>% lapply (hw.test, B = 1000))

#just extract p-values(repeat for each filter set)
(snps_hwe_pop_mat <- sapply(snps_hwe_pop, "[", i = TRUE, j = 3))
write.csv(x=snps_hwe_pop_mat, file = "snps_hwe_thinned_pop.csv")

################################################################################################

#diveRsity
#basic statistics (allelic richness etc using either rarefaction or bootstrapping)
#fis: inbreeding coefficient
#ar: allelic richness

#"snps_hwe_thinned_pop.txt" strict filter set with loci out of HWE removed, then converted to genepop 
#format in PGDspider 

basicStats(infile="snps_hwe_thinned_pop.txt", outfile = "strict_diveRsity_basicstats", 
           fis_ci = TRUE, ar_ci = TRUE, fis_boots = 1000, ar_boots = 1000,
           mc_reps = 10, rarefaction = TRUE, ar_alpha = 0.05,fis_alpha = 0.05)

divBasic(infile="snps_hwe_thinned_pop.txt", outfile = "strict_divbasic.txt", gp = 3, bootstraps = NULL,HWEexact = FALSE, mcRep = 2000)

#differentiation
fastDivPart(infile="snps_hwe_thinned_pop.txt", outfile = "strict_diveRsity_differentiation", pairwise = TRUE, fst = TRUE, bs_locus = FALSE,
            bs_pairwise = TRUE,boots = 100, plot = FALSE,para = FALSE)

#################################################################
#read in stats

#anova to compare pop genetic parameters
compare=read.csv("anova_mod.csv")

#friedman test instead:
f=friedman.test(Ar~region_time|stat, data=compare)
wilcox_test(Ar ~ region_time, paired = TRUE, p.adjust.method = "bonferroni", data=compare)
################

#observed heterozygosity
m1 = aov(obs_het~region*time+stat, data=compare)
summary(m1)
ho=TukeyHSD(m1)
ho

#expected heterozygosity
m2 = aov(exp_het~region*time+stat, data=compare)
summary(m2)
he=TukeyHSD(m2)
he

#################################################################
#DACP
#repeat for each filter set 

#NO a priori grouping 

#identify number of groups with k means 
grp1=find.clusters(rost_vcf_mod, max.n.clust=5)


#biggest increase after n=2
#look at group assignment 
table(pop(rost_vcf_mod), grp1$grp)

#cross-validation to determine number of principal components to use
#k-means: use as many PCs as need, DAPC: want to minimize PCs, overfitting bad 
set.seed(999)

rost1 <- xvalDapc(tab(rost_vcf_mod, NA.method = "mean"),  grp1$grp)

#check results 
rost1[-1]

#n.da is number of populations - 1 

dapc1 <- dapc(rost_gen02, var.contrib = TRUE, n.pca=25, n.da=1, grp1$grp)

#print contents of the object
print.dapc(dapc1)

#summary/useful info 
summary.dapc(dapc1)

scatter(dapc1,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

compoplot(dapc1, col=myCol,lab="", ncol=2)

loadingplot(dapc1$var.contr, threshold=quantile(dapc1$var.contr,0.75))

###############################################################
#DAPC: 

#a priori grouping 

#cross-validation to determine number of principal components to use
#k-means: use as many PCs as need, DAPC: want to minimize PCs, overfitting bad 
set.seed(999)
rostx <- xvalDapc(tab(rost_vcf_mod, NA.method = "mean"), pop(rost_vcf_mod))
#check results 
rostx[-1]

#n.da is number of populations - 1 
dapc2 <- dapc(rost_vcf_mod, var.contrib = TRUE, n.pca=25, n.da=3, pop(rost_vcf_mod))

#print contents of the object
print.dapc(dapc2)
#summary/useful info 
summary.dapc(dapc2)
#predict individual assignment 
predict.dapc(dapc2)

scatter(dapc2,scree.da=FALSE, scree.pca=TRUE, bg="white", posi.pca="topleft", 
        legend=TRUE, col=myCol, clab=0, cstar=0, cex=2, pch=c(15, 16, 17, 18), solid=0.8)

compoplot(dapc2, col=myCol,lab="", ncol=2)

loadingplot(dapc2$var.contr, threshold=quantile(dapc2$var.contr,0.75))

#plot of population assignment predict.dapc(dapc2)
assignplot(dapc2, subset=1:76)
