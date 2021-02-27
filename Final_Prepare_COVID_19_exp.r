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
Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_WBC_01_21/"

 
load("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/SNPS_37_38.RData")

for (n in c("HGI_round_4_A2","HGI_round_4_B2","HGI_round_5_A2","HGI_round_5_B2","HGI_round_5_C2")) {
  rm(exp_dat)
  rm(exp_dat1)
  exp_dat_file=paste(Pathway_SNP, n, "_38.txt", sep="")
  exp_dat <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
  names(exp_dat)[1] <- "CHR"
  exp_dat$location_38=paste(exp_dat$CHR, exp_dat$POS , sep=":")
  exp_dat<-exp_dat %>% left_join(SNPS_37_38, by= "location_38")
  exp_dat1=exp_dat[exp_dat$all_inv_var_meta_p<5e-8, ]
  Outputfile=paste(Pathway_SNP, n, "_as_exposure.txt", sep="")
  write.table(exp_dat, file= Outputfile, append = TRUE, row.names = FALSE,col.names = T, quote = FALSE, sep='\t')
  Outputfile1=paste(Pathway_SNP, n, "_as_exposure_exp_5e-8.txt", sep="")
  write.table(exp_dat1, file= Outputfile1, append = TRUE, row.names = FALSE,col.names = T, quote = FALSE, sep='\t')
}
           


for (n in c("HGI_round_4_A2","HGI_round_4_B2","HGI_round_5_A2","HGI_round_5_B2","HGI_round_5_C2")) {
  outcomefile=paste(Pathway_SNP,  n ,"_as_exposure_exp_5e-8.txt", sep="")
  ####### Change csv file
  exp_dat2 <- read_exposure_data(
    filename = outcomefile,
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "all_inv_var_meta_beta",
    se_col = "all_inv_var_meta_sebeta",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "all_meta_AF",
    samplesize_col = "all_meta_sample_N",
    pval_col = "all_inv_var_meta_p"
    #id_col = Trait$V1[e],
    #phenotype_col = Trait$V1[e],
    #units_col = "Units",
    #gene_col = "Gene",
  )
  exp_dat2$id.exposure = n
  exp_dat2$exposure = n
  
  exp_dat3=clump_data(
    exp_dat2,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR"
  )
  colnames(exp_dat3)[colnames(exp_dat3)=="SNP"] <- "variant_id"
  print(exp_dat3$variant_id)
  COVID <- read.csv(outcomefile,header=T, as.is=T,sep = "\t")
  exp_dat4 <- exp_dat3 %>% inner_join(COVID, by= "variant_id")
  Outputfile=paste(Pathway_SNP, n ,"_SNP.txt", sep="")
  write.table(exp_dat4, file= Outputfile, col.names = TRUE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
}


