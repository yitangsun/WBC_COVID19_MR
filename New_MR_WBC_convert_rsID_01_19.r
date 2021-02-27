# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")

library(TwoSampleMR)
library(dplyr)
library(googlesheets)
library(vcfR)
library(stringr)
library(mr.raps)
'%ni%' <- Negate('%in%')
# ieugwasr::api_status()
# $`API version`
# #3.6.7
# $`Total associations`
# [1] 126335269652
# $`Total complete datasets`
# [1] 34670
# $`Total public datasets`
# [1] 34513
library(MRInstruments)
library(MVMR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library("SNPlocs.Hsapiens.dbSNP151.GRCh38")


ao <- available_outcomes()
# data(gwas_catalog)
# head(gwas_catalog)

Pathway_SNP="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/"
Pathway_geno="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/"
Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_all_10_9/F_stat/"

Trait <- read.csv("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/WBC_Trait.txt",header=F, as.is=T,sep = "\t")
Trait$V4=as.numeric(gsub(",", "", Trait$V4))

load("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/SNPS_37.RData")
#load("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/SNPS_38.RData")
#load("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/SNPS_37_38.RData")

len_exp_file=length(Trait$V1)
for (e in 1:len_exp_file) {
  if (Trait$V2[e]==1) {
    exp_dat_file=paste(Pathway_SNP, Trait$V1[e], "_buildGRCh37.tsv", sep="")
    exp_dat_temtem <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
    exp_dat= exp_dat_temtem[exp_dat_temtem$p_value<5e-8, ]
    Outputfile=paste(Pathway_SNP, Trait$V1[e], "_5e-8", sep="")
    write.table(exp_dat, file= Outputfile, col.names = TRUE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  } else if (Trait$V2[e]==2) {
    exp_dat_file=paste(Pathway_SNP, Trait$V1[e], "_buildGRCh37.tsv", sep="")
    exp_dat_temtem <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
    exp_dat= exp_dat_temtem[exp_dat_temtem$p_value<5e-8, ]
    exp_dat$location_37=paste(exp_dat$chromosome,":",exp_dat$base_pair_location, sep="")
    exp_dat <- exp_dat %>% left_join(SNPS_37, by= "location_37")
    Outputfile=paste(Pathway_SNP, Trait$V1[e], "_5e-8", sep="")
    write.table(exp_dat, file= Outputfile, col.names = TRUE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  } else if (Trait$V2[e]==3) {
    exp_dat_file=paste(Pathway_SNP,"27863252-", Trait$V1[e], ".tsv", sep="")
    exp_dat_temtem <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
    exp_dat= exp_dat_temtem[exp_dat_temtem$p_value<5e-8, ]
    Outputfile=paste(Pathway_SNP, Trait$V1[e], "_5e-8", sep="")
    write.table(exp_dat, file= Outputfile, col.names = TRUE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  }
}

#######Manually change rsID in Trait$V2[e]==2 _5e-8       
