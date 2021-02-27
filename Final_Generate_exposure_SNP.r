
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

Pathway_SNP="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/"
Pathway_geno="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/"
Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_WBC_01_28_5e-8/"

#load("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/SNPS_37_38.RData")

Trait <- read.csv("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/WBC_Trait.txt",header=F, as.is=T,sep = "\t")
Trait$V4=as.numeric(gsub(",", "", Trait$V4))

len_exp_file=length(Trait$V1)
for (e in 1:len_exp_file) {
  if (Trait$V2[e]==0) {
    exp_dat <- extract_instruments(Trait$V1[e])
    exp_dat$pos.exposure_end=exp_dat$pos.exposure+1
    exp_dat$chr.exposure_new=paste("chr",exp_dat$chr.exposure, sep="") 
    exp_dat<-exp_dat%>%select(chr.exposure_new, pos.exposure, pos.exposure_end, SNP, chr.exposure)
    write.table(exp_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt", append = TRUE, row.names = FALSE,col.names = FALSE, quote = FALSE, sep=' ')
  } else if (Trait$V2[e]==10) {
    exp_dat <- extract_instruments(Trait$V1[e])
    exp_dat$pos.exposure_end=exp_dat$pos.exposure+1
    exp_dat$chr.exposure_new=paste("chr",exp_dat$chr.exposure, sep="") 
    exp_dat<-exp_dat%>%select(chr.exposure_new, pos.exposure, pos.exposure_end, SNP, chr.exposure)
    write.table(exp_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt", append = TRUE, row.names = FALSE,col.names = FALSE, quote = FALSE, sep=' ')
  } else if (Trait$V2[e]==30) {
    exp_dat <- extract_instruments(Trait$V1[e])
    exp_dat$pos.exposure_end=exp_dat$pos.exposure+1
    exp_dat$chr.exposure_new=paste("chr",exp_dat$chr.exposure, sep="") 
    exp_dat<-exp_dat%>%select(chr.exposure_new, pos.exposure, pos.exposure_end, SNP, chr.exposure)
    write.table(exp_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt", append = TRUE, row.names = FALSE,col.names = FALSE, quote = FALSE, sep=' ')
  } else if (Trait$V2[e]==1) {
    exp_dat_file=paste(Pathway_SNP, Trait$V1[e], "_5e-8", sep="")
    exp_dat <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
    colnames(exp_dat)[colnames(exp_dat)=="variant_id"] <- "SNP"
    colnames(exp_dat)[colnames(exp_dat)=="base_pair_location"] <- "pos.exposure"
    colnames(exp_dat)[colnames(exp_dat)=="chromosome"] <- "chr.exposure"
    exp_dat$pos.exposure_end=exp_dat$pos.exposure+1
    exp_dat$chr.exposure_new=paste("chr",exp_dat$chr.exposure, sep="") 
    exp_dat<-exp_dat%>%select(chr.exposure_new, pos.exposure, pos.exposure_end, SNP, chr.exposure)
    write.table(exp_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt", append = TRUE, row.names = FALSE,col.names = FALSE, quote = FALSE, sep=' ')
  } else if (Trait$V2[e]==2) {
    exp_dat_file=paste(Pathway_SNP, Trait$V1[e], "_5e-8", sep="")
    exp_dat <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
    colnames(exp_dat)[colnames(exp_dat)=="variant_id"] <- "SNP"
    colnames(exp_dat)[colnames(exp_dat)=="base_pair_location"] <- "pos.exposure"
    colnames(exp_dat)[colnames(exp_dat)=="chromosome"] <- "chr.exposure"
    exp_dat$pos.exposure_end=exp_dat$pos.exposure+1
    exp_dat$chr.exposure_new=paste("chr",exp_dat$chr.exposure, sep="") 
    exp_dat<-exp_dat%>%select(chr.exposure_new, pos.exposure, pos.exposure_end, SNP, chr.exposure)
    write.table(exp_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt", append = TRUE, row.names = FALSE,col.names = FALSE, quote = FALSE, sep=' ')
  } else if (Trait$V2[e]==3) {
    exp_dat_file=paste(Pathway_SNP, Trait$V1[e], "_5e-8", sep="")
    exp_dat <- read.csv(exp_dat_file,header=T, as.is=T,sep = "\t")
    exp_dat$location_38=paste(exp_dat$hm_chrom, exp_dat$hm_pos , sep=":")
    exp_dat<-exp_dat %>% left_join(SNPS_37_38, by= "location_38")
    colnames(exp_dat)[colnames(exp_dat)=="hm_rsid"] <- "SNP"
    colnames(exp_dat)[colnames(exp_dat)=="pos_37"] <- "pos.exposure"
    colnames(exp_dat)[colnames(exp_dat)=="chr_37"] <- "chr.exposure"
    exp_dat$pos.exposure_end=exp_dat$pos.exposure+1
    exp_dat$chr.exposure_new=paste("chr",exp_dat$chr.exposure, sep="") 
    exp_dat<-exp_dat%>%select(chr.exposure_new, pos.exposure, pos.exposure_end, SNP, chr.exposure)
    write.table(exp_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt", append = TRUE, row.names = FALSE,col.names = FALSE, quote = FALSE, sep=' ')
  }
}


######### Generate outcome dataset using Exposure_SNP_WBC.txt

tem_dat=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/HGI_round_4_A2.txt',header = T, as.is=T, sep='\t')
exp_dat_SNP=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt',header = F, as.is=T, sep=' ')
exp_dat_SNP=exp_dat_SNP[!duplicated(exp_dat_SNP$V4), ]

# exp_dat_SNP$V1=as.character.factor(exp_dat_SNP$V1)
# exp_dat_SNP$V4=as.character.factor(exp_dat_SNP$V4)
# exp_dat_SNP$V5=as.character.factor(exp_dat_SNP$V5)
exp_dat_SNP=exp_dat_SNP[exp_dat_SNP$V5 !="X",]
exp_dat_SNP$V5=as.integer(exp_dat_SNP$V5)
colnames(tem_dat)[1] <- "CHR"
tem_dat$SNP=as.integer(tem_dat$SNP)

tem_dat$SNP_tem=paste(tem_dat$CHR, tem_dat$POS , sep="_")
exp_dat_SNP$SNP_tem=paste(exp_dat_SNP$V5, exp_dat_SNP$V2 , sep="_")

tem_dat<-exp_dat_SNP %>% left_join(tem_dat, by= "SNP_tem")
write.table(tem_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/HGI_round_4_A2.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep='\t')


tem_dat=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/HGI_round_4_B2.txt',header = T, as.is=T, sep='\t')
exp_dat_SNP=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt',header = F, as.is=T, sep=' ')
exp_dat_SNP=exp_dat_SNP[!duplicated(exp_dat_SNP$V4), ]

# exp_dat_SNP$V1=as.character.factor(exp_dat_SNP$V1)
# exp_dat_SNP$V4=as.character.factor(exp_dat_SNP$V4)
# exp_dat_SNP$V5=as.character.factor(exp_dat_SNP$V5)
exp_dat_SNP=exp_dat_SNP[exp_dat_SNP$V5 !="X",]
exp_dat_SNP$V5=as.integer(exp_dat_SNP$V5)
colnames(tem_dat)[1] <- "CHR"
tem_dat$SNP=as.integer(tem_dat$SNP)

tem_dat$SNP_tem=paste(tem_dat$CHR, tem_dat$POS , sep="_")
exp_dat_SNP$SNP_tem=paste(exp_dat_SNP$V5, exp_dat_SNP$V2 , sep="_")

tem_dat<-exp_dat_SNP %>% left_join(tem_dat, by= "SNP_tem")
write.table(tem_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/HGI_round_4_B2.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep='\t')



tem_dat=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/HGI_round_5_A2.txt',header = T, as.is=T, sep='\t')
exp_dat_SNP=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt',header = F, as.is=T, sep=' ')
exp_dat_SNP=exp_dat_SNP[!duplicated(exp_dat_SNP$V4), ]

# exp_dat_SNP$V1=as.character.factor(exp_dat_SNP$V1)
# exp_dat_SNP$V4=as.character.factor(exp_dat_SNP$V4)
# exp_dat_SNP$V5=as.character.factor(exp_dat_SNP$V5)
exp_dat_SNP=exp_dat_SNP[exp_dat_SNP$V5 !="X",]
exp_dat_SNP$V5=as.integer(exp_dat_SNP$V5)
colnames(tem_dat)[1] <- "CHR"
tem_dat$SNP=as.integer(tem_dat$SNP)

tem_dat$SNP_tem=paste(tem_dat$CHR, tem_dat$POS , sep="_")
exp_dat_SNP$SNP_tem=paste(exp_dat_SNP$V5, exp_dat_SNP$V2 , sep="_")

tem_dat<-exp_dat_SNP %>% left_join(tem_dat, by= "SNP_tem")
write.table(tem_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/HGI_round_5_A2.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep='\t')



tem_dat=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/HGI_round_5_B2.txt',header = T, as.is=T, sep='\t')
exp_dat_SNP=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt',header = F, as.is=T, sep=' ')
exp_dat_SNP=exp_dat_SNP[!duplicated(exp_dat_SNP$V4), ]

# exp_dat_SNP$V1=as.character.factor(exp_dat_SNP$V1)
# exp_dat_SNP$V4=as.character.factor(exp_dat_SNP$V4)
# exp_dat_SNP$V5=as.character.factor(exp_dat_SNP$V5)
exp_dat_SNP=exp_dat_SNP[exp_dat_SNP$V5 !="X",]
exp_dat_SNP$V5=as.integer(exp_dat_SNP$V5)
colnames(tem_dat)[1] <- "CHR"
tem_dat$SNP=as.integer(tem_dat$SNP)

tem_dat$SNP_tem=paste(tem_dat$CHR, tem_dat$POS , sep="_")
exp_dat_SNP$SNP_tem=paste(exp_dat_SNP$V5, exp_dat_SNP$V2 , sep="_")

tem_dat<-exp_dat_SNP %>% left_join(tem_dat, by= "SNP_tem")
write.table(tem_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/HGI_round_5_B2.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep='\t')


tem_dat=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/HGI_round_5_C2.txt',header = T, as.is=T, sep='\t')
exp_dat_SNP=read.csv('/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/Exposure_SNP_WBC.txt',header = F, as.is=T, sep=' ')
exp_dat_SNP=exp_dat_SNP[!duplicated(exp_dat_SNP$V4), ]

# exp_dat_SNP$V1=as.character.factor(exp_dat_SNP$V1)
# exp_dat_SNP$V4=as.character.factor(exp_dat_SNP$V4)
# exp_dat_SNP$V5=as.character.factor(exp_dat_SNP$V5)
exp_dat_SNP=exp_dat_SNP[exp_dat_SNP$V5 !="X",]
exp_dat_SNP$V5=as.integer(exp_dat_SNP$V5)
colnames(tem_dat)[1] <- "CHR"
tem_dat$SNP=as.integer(tem_dat$SNP)

tem_dat$SNP_tem=paste(tem_dat$CHR, tem_dat$POS , sep="_")
exp_dat_SNP$SNP_tem=paste(exp_dat_SNP$V5, exp_dat_SNP$V2 , sep="_")

tem_dat<-exp_dat_SNP %>% left_join(tem_dat, by= "SNP_tem")
write.table(tem_dat, file= "/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/Genotype/HGI_round_5_C2.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep='\t')



   
