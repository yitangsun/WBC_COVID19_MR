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

Outcome_list=c("HGI_round_4_A2","HGI_round_4_B2","HGI_round_5_A2","HGI_round_5_B2","HGI_round_5_C2")


########Function
My_MR <- function(exp_dat,outcome_dat) {
  rm(dat)
  try(dat<- exp_dat %>% inner_join(outcome_dat, by= "SNP"), silent=TRUE)
  dat=dat[is.na(dat$effect_allele.exposure)==FALSE ,]
  dat=dat[is.na(dat$effect_allele.outcome)==FALSE ,]
  dat=dat[is.na(dat$other_allele.exposure)==FALSE ,]
  dat=dat[is.na(dat$other_allele.outcome)==FALSE ,]
  if (exists("dat")==TRUE && length(dat$SNP)>0) {
    for (y in 1:length(dat$SNP)) {
      if ( dat$effect_allele.exposure[y]==dat$effect_allele.outcome[y] && dat$other_allele.exposure[y]==dat$other_allele.outcome[y] ) {
        dat$betaYG[y] =dat$beta.outcome[y]
      } else if (dat$effect_allele.exposure[y]==dat$other_allele.outcome[y] && dat$other_allele.exposure[y]==dat$effect_allele.outcome[y]) {
        dat$betaYG[y] =-1*dat$beta.outcome[y]
      } else{
        dat$betaYG[y]= NA
      }
    }
    dat$mr_keep=TRUE
    colnames(dat)[colnames(dat)=="beta.outcome"] <- "beta.o_old"
    colnames(dat)[colnames(dat)=="betaYG"] <- "beta.outcome"
    
    dat=dat[is.na(dat$beta.outcome)==FALSE ,]
    
    dat$len_SNP=length(dat$SNP)
    dat$R_Square = 2*(dat$beta.exposure^2)*dat$eaf.exposure*(1-dat$eaf.exposure)
    dat$Total_R_Square = sum(dat$R_Square)
    dat$F_stat=dat$Total_R_Square*(dat$len_SNP)*(dat$samplesize.exposure-dat$len_SNP-1)/(1-dat$Total_R_Square)
    dat$F_stat_sim=dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP)
    F_stat=data.frame(dat$id.exposure[1],dat$len_SNP[1],dat$samplesize.exposure[1],dat$Total_R_Square[1],dat$F_stat[1],dat$F_stat_sim[1])
    colnames(F_stat) <- c("id.exposure","len_SNP","samplesize.exposure","Total_R_Square","F_stat","F_stat_sim")
    
    try(dat<-dat%>%select(beta.exposure, se.exposure, beta.outcome,se.outcome, mr_keep, id.exposure, id.outcome, 
                          exposure, outcome, SNP,pval.exposure, pval.outcome, samplesize.exposure, samplesize.outcome,
                          len_SNP, R_Square, Total_R_Square, F_stat, F_stat_sim), silent=TRUE)
    
    # dat$af=0.07
    # dat$ncase=1610
    # dat$ncontrol=2180
    # dat$prevalence=0.05
    # 
    # dat$b=get_r_from_lor (lor=dat$beta.outcome,af=dat$af,ncase=dat$ncase,ncontrol=dat$ncontrol,prevalence=dat$prevalence)
    # 
    # snp_r2.exposure=mr_steiger(p_exp=dat$pval.exposure, p_out=dat$pval.outcome, n_exp=dat$samplesize.exposure, n_out=dat$samplesize_col, r_exp=dat$a, r_out=dat$b)$r2_exp
    # snp_r2.outcome=mr_steiger(p_exp=dat$pval.exposure, p_out=dat$pval.outcome, n_exp=dat$samplesize.exposure, n_out=dat$samplesize_col, r_exp=dat$a, r_out=dat$b)$r2_out
    # Correct_causal_direction=mr_steiger(p_exp=dat$pval.exposure, p_out=dat$pval.outcome, n_exp=dat$samplesize.exposure, n_out=dat$samplesize_col, r_exp=dat$a, r_out=dat$b)$correct_causal_direction
    # steiger_pval=mr_steiger(p_exp=dat$pval.exposure, p_out=dat$pval.outcome, n_exp=dat$samplesize.exposure, n_out=dat$samplesize_col, r_exp=dat$a, r_out=dat$b)$steiger_test
    # Steiger=data.frame(dat$id.exposure[1], snp_r2.exposure, snp_r2.outcome, Correct_causal_direction, steiger_pval)
    # colnames(Steiger)[1] <- "id.exposure"
    
    if (length(which(dat$mr_keep=='TRUE'))>3) {
      PressoObject=run_mr_presso(dat, NbDistribution = length(dat$SNP)/0.05, SignifThreshold = 0.05)
      if (length(PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]) >0 && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]!=c("No significant outliers") && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]!=c("All SNPs considered as outliers")) {
        # OutliersIndices=PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]
        # dat=dat[-OutliersIndices, ]
        if (length(which(dat$mr_keep=='TRUE'))>2) {
          pleiotropy=mr_pleiotropy_test(dat)
          if (pleiotropy$pval<0.05) {
            Heterogeneity=mr_heterogeneity(dat)
            res <- mr(dat)
            final_res = res[res$method=='MR Egger',]
          } else {
            Heterogeneity=mr_heterogeneity(dat)
            if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval<0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval<0.05) {
              if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval < Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval) {
                res <- mr(dat)
                final_res = res[res$method=='MR Egger',]
              } else {
                res <- mr(dat)
                final_res = res[res$method=='Inverse variance weighted',]
              }
            } else if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval>=0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval<0.05) {
              res <- mr(dat)
              final_res = res[res$method=='Inverse variance weighted',]
            } else if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval<0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval>=0.05) {
              res <- mr(dat)
              final_res = res[res$method=='MR Egger',]
            } else{
              res <- mr(dat)
              final_res = res[res$method=='Inverse variance weighted',]
            }
          }
        } else if (length(which(dat$mr_keep=='TRUE'))==1) {
          res <- mr(dat)
          final_res = res[res$method=='Wald ratio',]
        } else if (length(which(dat$mr_keep=='TRUE'))==2) {
          res <- mr(dat)
          final_res = res[res$method=='Inverse variance weighted',]
        }
      } else if (length(PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]) >0 && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]!=c("No significant outliers") && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]==c("All SNPs considered as outliers")) {
        Heterogeneity=mr_heterogeneity(dat)
        res <- mr(dat)
        final_res = res[res$method=='MR Egger',]
      } else {
        if (length(which(dat$mr_keep=='TRUE'))>2) {
          #Horizontal pleiotropy
          pleiotropy=mr_pleiotropy_test(dat)
          if (pleiotropy$pval<0.05) {
            Heterogeneity=mr_heterogeneity(dat)
            res <- mr(dat)
            final_res = res[res$method=='MR Egger',]
          } else {
            # Heterogeneity statistics
            Heterogeneity=mr_heterogeneity(dat)
            if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval<0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval<0.05) {
              if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval < Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval) {
                res <- mr(dat)
                final_res = res[res$method=='MR Egger',]
              } else {
                res <- mr(dat)
                final_res = res[res$method=='Inverse variance weighted',]
              }
            } else if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval>=0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval<0.05) {
              res <- mr(dat)
              final_res = res[res$method=='Inverse variance weighted',]
            } else if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval<0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval>=0.05) {
              res <- mr(dat)
              final_res = res[res$method=='MR Egger',]
            } else{
              res <- mr(dat)
              final_res = res[res$method=='Inverse variance weighted',]
            }
          }
        } else if (length(which(dat$mr_keep=='TRUE'))==1) {
          res <- mr(dat)
          final_res = res[res$method=='Wald ratio',]
        } else if (length(which(dat$mr_keep=='TRUE'))==2) {
          res <- mr(dat)
          final_res = res[res$method=='Inverse variance weighted',]
        }
      } 
    } else {
      if (length(which(dat$mr_keep=='TRUE'))>2) {
        #Horizontal pleiotropy
        pleiotropy=mr_pleiotropy_test(dat)
        if (pleiotropy$pval<0.05) {
          Heterogeneity=mr_heterogeneity(dat)
          res <- mr(dat)
          final_res = res[res$method=='MR Egger',]
        } else {
          # Heterogeneity statistics
          Heterogeneity=mr_heterogeneity(dat)
          if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval<0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval<0.05) {
            if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval < Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval) {
              res <- mr(dat)
              final_res = res[res$method=='MR Egger',]
            } else {
              res <- mr(dat)
              final_res = res[res$method=='Inverse variance weighted',]
            }
          } else if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval>=0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval<0.05) {
            res <- mr(dat)
            final_res = res[res$method=='Inverse variance weighted',]
          } else if (Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]$Q_pval<0.05 && Heterogeneity[Heterogeneity$method=='MR Egger',]$Q_pval>=0.05) {
            res <- mr(dat)
            final_res = res[res$method=='MR Egger',]
          } else{
            res <- mr(dat)
            final_res = res[res$method=='Inverse variance weighted',]
          }
        }
      } else if (length(which(dat$mr_keep=='TRUE'))==1) {
        res <- mr(dat)
        final_res = res[res$method=='Wald ratio',]
      } else if (length(which(dat$mr_keep=='TRUE'))==2) {
        res <- mr(dat)
        final_res = res[res$method=='Inverse variance weighted',]
      }
    }
  }
  
  rm(dat)
  try(dat<- exp_dat %>% inner_join(outcome_dat, by= "SNP"), silent=TRUE)
  dat=dat[is.na(dat$effect_allele.exposure)==FALSE ,]
  dat=dat[is.na(dat$effect_allele.outcome)==FALSE ,]
  dat=dat[is.na(dat$other_allele.exposure)==FALSE ,]
  dat=dat[is.na(dat$other_allele.outcome)==FALSE ,]
  if (exists("dat")==TRUE && length(dat$SNP)>0) {
    for (y in 1:length(dat$SNP)) {
      if ( dat$effect_allele.exposure[y]==dat$effect_allele.outcome[y] && dat$other_allele.exposure[y]==dat$other_allele.outcome[y] ) {
        dat$betaYG[y] =dat$beta.outcome[y]
      } else if (dat$effect_allele.exposure[y]==dat$other_allele.outcome[y] && dat$other_allele.exposure[y]==dat$effect_allele.outcome[y]) {
        dat$betaYG[y] =-1*dat$beta.outcome[y]
      } else{
        dat$betaYG[y]= NA
      }
    }
    dat$mr_keep=TRUE
    colnames(dat)[colnames(dat)=="beta.outcome"] <- "beta.o_old"
    colnames(dat)[colnames(dat)=="betaYG"] <- "beta.outcome"
    try(dat<-dat%>%select(beta.exposure, se.exposure, beta.outcome,se.outcome, mr_keep, id.exposure, id.outcome, 
                          exposure, outcome, SNP,pval.exposure, pval.outcome, samplesize.exposure, samplesize.outcome), silent=TRUE)
    
    dat=dat[is.na(dat$beta.outcome)==FALSE ,]
    
    if (length(which(dat$mr_keep=='TRUE'))>3) {
      PressoObject=run_mr_presso(dat, NbDistribution = length(dat$SNP)/0.05, SignifThreshold = 0.05)
      if (length(PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]) >0 && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]!=c("No significant outliers") && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]!=c("All SNPs considered as outliers")) {
        # OutliersIndices=PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]
        # dat=dat[-OutliersIndices, ]
        if (length(which(dat$mr_keep=='TRUE'))>2) {
          #Horizontal pleiotropy
          pleiotropy=mr_pleiotropy_test(dat)
          res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw","mr_egger_regression", "mr_two_sample_ml", "mr_simple_median", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
          IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
          MR_Egger = res[res$method=='MR Egger',]
          IVW = res[res$method=='Inverse variance weighted',]
          ml = res[res$method=='Maximum likelihood',]
          W_Med = res[res$method=='Weighted median',]
          W_Mod = res[res$method=='Weighted mode',]
        } else if (length(which(dat$mr_keep=='TRUE'))==1) {
          res <- mr(dat, method_list=c("mr_wald_ratio"))
          Wald = res[res$method=='Wald ratio',]
        } else if (length(which(dat$mr_keep=='TRUE'))==2) {
          res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw", "mr_two_sample_ml"))
          IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
          IVW = res[res$method=='Inverse variance weighted',]
          ml = res[res$method=='Maximum likelihood',]
        }
      } else if (length(PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]) >0 && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]!=c("No significant outliers") && PressoObject[[1]][['MR-PRESSO results']][['Distortion Test']][['Outliers Indices']]==c("All SNPs considered as outliers")) {
        if (length(which(dat$mr_keep=='TRUE'))>2) {
          #Horizontal pleiotropy
          pleiotropy=mr_pleiotropy_test(dat)
          res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw","mr_egger_regression", "mr_two_sample_ml", "mr_simple_median", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
          IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
          MR_Egger = res[res$method=='MR Egger',]
          IVW = res[res$method=='Inverse variance weighted',]
          ml = res[res$method=='Maximum likelihood',]
          W_Med = res[res$method=='Weighted median',]
          W_Mod = res[res$method=='Weighted mode',]
        } else if (length(which(dat$mr_keep=='TRUE'))==1) {
          res <- mr(dat, method_list=c("mr_wald_ratio"))
          Wald = res[res$method=='Wald ratio',]
        } else if (length(which(dat$mr_keep=='TRUE'))==2) {
          res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw", "mr_two_sample_ml"))
          IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
          IVW = res[res$method=='Inverse variance weighted',]
          ml = res[res$method=='Maximum likelihood',]
        }
      } else {
        if (length(which(dat$mr_keep=='TRUE'))>2) {
          #Horizontal pleiotropy
          pleiotropy=mr_pleiotropy_test(dat)
          res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw","mr_egger_regression", "mr_two_sample_ml", "mr_simple_median", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
          IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
          MR_Egger = res[res$method=='MR Egger',]
          IVW = res[res$method=='Inverse variance weighted',]
          ml = res[res$method=='Maximum likelihood',]
          W_Med = res[res$method=='Weighted median',]
          W_Mod = res[res$method=='Weighted mode',]
        } else if (length(which(dat$mr_keep=='TRUE'))==1) {
          res <- mr(dat, method_list=c("mr_wald_ratio"))
          Wald = res[res$method=='Wald ratio',]
        } else if (length(which(dat$mr_keep=='TRUE'))==2) {
          res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw", "mr_two_sample_ml"))
          IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
          IVW = res[res$method=='Inverse variance weighted',]
          ml = res[res$method=='Maximum likelihood',]
        }
      }
    } else {
      if (length(which(dat$mr_keep=='TRUE'))>2) {
        #Horizontal pleiotropy
        pleiotropy=mr_pleiotropy_test(dat)
        res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw","mr_egger_regression", "mr_two_sample_ml", "mr_simple_median", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))
        MR_Egger = res[res$method=='MR Egger',]
        IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
        IVW = res[res$method=='Inverse variance weighted',]
        ml = res[res$method=='Maximum likelihood',]
        W_Med = res[res$method=='Weighted median',]
        W_Mod = res[res$method=='Weighted mode',]
      } else if (length(which(dat$mr_keep=='TRUE'))==1) {
        res <- mr(dat, method_list=c("mr_wald_ratio"))
        Wald = res[res$method=='Wald ratio',]
      } else if (length(which(dat$mr_keep=='TRUE'))==2) {
        res <- mr(dat, method_list=c("mr_ivw_mre","mr_ivw", "mr_two_sample_ml"))
        IVW_MRE = res[res$method=='Inverse variance weighted (multiplicative random effects)',]
        IVW = res[res$method=='Inverse variance weighted',]
        ml = res[res$method=='Maximum likelihood',]
      }
    }
  }
  
  if (exists("final_res")==TRUE ) {
    final_res$Test=n
    names(final_res)[names(final_res) == "method"] <- "Methods"
    names(final_res)[names(final_res) == "nsnp"] <- "nsnps"
    names(final_res)[names(final_res) == "b"] <- "Beta"
    names(final_res)[names(final_res) == "se"] <- "SE"
    names(final_res)[names(final_res) == "exposure"] <- "trait"
    
    if (exists("F_stat")==TRUE){
      final_res<- F_stat %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$len_SNP=NA
      final_res$samplesize.exposure=NA
      final_res$Total_R_Square=NA
      final_res$F_stat=NA
      final_res$F_stat_sim=NA
    }
    
    # if (exists("Steiger")==TRUE){
    #   final_res<- Steiger %>% right_join(final_res, by= "id.exposure")
    # } else {
    #   final_res$snp_r2.exposure=NA
    #   final_res$snp_r2.outcome=NA
    #   final_res$Correct_causal_direction=NA
    #   final_res$steiger_pval=NA
    # }
    
    if (exists("Heterogeneity")==TRUE){
      Het_IVW=Heterogeneity[Heterogeneity$method=='Inverse variance weighted',]
      names(Het_IVW)[names(Het_IVW) == "Q_pval"] <- "Het_IVW_pval"
      final_res<- Het_IVW %>% right_join(final_res, by= "id.exposure")
      Het_MR_Egger=Heterogeneity[Heterogeneity$method=='MR Egger',]
      names(Het_MR_Egger)[names(Het_MR_Egger) == "Q_pval"] <- "Het_Egger_pval"
      final_res<- Het_MR_Egger %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$Het_IVW_pval=NA
      final_res$Het_Egger_pval=NA
    }
    
    if (exists("pleiotropy")==TRUE){
      names(pleiotropy)[names(pleiotropy) == "egger_intercept"] <- "Egger_intercept"
      names(pleiotropy)[names(pleiotropy) == "pval"] <- "pval_intercept"
      final_res<- pleiotropy %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$Egger_intercept=NA
      final_res$pval_intercept=NA
    }
    
    if (exists("PressoObject")==TRUE){
      PRESSO = PressoObject[[1]]$`Main MR results`
      PRESSO_raw = PRESSO[PRESSO$"MR Analysis"=='Raw',]
      PRESSO_corrected = PRESSO[PRESSO$"MR Analysis"=='Outlier-corrected',]
      PRESSO_GLOBAL = PressoObject[[1]]$`MR-PRESSO results`$`Global Test`
      pval_PRESSO_Global = PRESSO_GLOBAL$Pvalue
      names(PRESSO_raw)[names(PRESSO_raw) == "Causal Estimate"] <- "b_PRESSO_raw"
      names(PRESSO_raw)[names(PRESSO_raw) == "Sd"] <- "se_PRESSO_raw"
      names(PRESSO_raw)[names(PRESSO_raw) == "P-value"] <- "pval_PRESSO_raw"
      names(PRESSO_corrected)[names(PRESSO_corrected) == "Causal Estimate"] <- "b_PRESSO_corrected"
      names(PRESSO_corrected)[names(PRESSO_corrected) == "Sd"] <- "se_PRESSO_corrected"
      names(PRESSO_corrected)[names(PRESSO_corrected) == "P-value"] <- "pval_PRESSO_corrected"
      PRESSO_raw<-PRESSO_raw%>%select(b_PRESSO_raw, se_PRESSO_raw, pval_PRESSO_raw)
      PRESSO_corrected<-PRESSO_corrected%>%select(b_PRESSO_corrected, se_PRESSO_corrected, pval_PRESSO_corrected)
      final_res<- cbind(PRESSO_corrected, final_res)
      final_res<- cbind(PRESSO_raw, final_res)
      final_res<- cbind(pval_PRESSO_Global, final_res)
    } else {
      final_res$b_PRESSO_raw=NA
      final_res$se_PRESSO_raw=NA
      final_res$pval_PRESSO_raw=NA
      final_res$b_PRESSO_corrected=NA
      final_res$se_PRESSO_corrected=NA
      final_res$pval_PRESSO_corrected=NA
      final_res$pval_PRESSO_Global=NA
    }
    
    if (exists("W_Mod")==TRUE){
      names(W_Mod)[names(W_Mod) == "b"] <- "b_W_Mod"
      names(W_Mod)[names(W_Mod) == "se"] <- "se_W_Mod"
      names(W_Mod)[names(W_Mod) == "pval"] <- "pval_W_Mod"
      final_res<- W_Mod %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_W_Mod=NA
      final_res$se_W_Mod=NA
      final_res$pval_W_Mod=NA
    }
    
    if (exists("W_Med")==TRUE){
      names(W_Med)[names(W_Med) == "b"] <- "b_W_Med"
      names(W_Med)[names(W_Med) == "se"] <- "se_W_Med"
      names(W_Med)[names(W_Med) == "pval"] <- "pval_W_Med"
      final_res<- W_Med %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_W_Med=NA
      final_res$se_W_Med=NA
      final_res$pval_W_Med=NA
    }
    
    if (exists("ml")==TRUE){
      names(ml)[names(ml) == "b"] <- "b_ml"
      names(ml)[names(ml) == "se"] <- "se_ml"
      names(ml)[names(ml) == "pval"] <- "pval_ml"
      final_res<- ml %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_ml=NA
      final_res$se_ml=NA
      final_res$pval_ml=NA
    }
    
    if (exists("Wald")==TRUE){
      names(Wald)[names(Wald) == "b"] <- "b_Wald"
      names(Wald)[names(Wald) == "se"] <- "se_Wald"
      names(Wald)[names(Wald) == "pval"] <- "pval_Wald"
      final_res<- Wald %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_Wald=NA
      final_res$se_Wald=NA
      final_res$pval_Wald=NA
    }
    
    if (exists("MR_Egger")==TRUE){
      names(MR_Egger)[names(MR_Egger) == "b"] <- "b_Egger"
      names(MR_Egger)[names(MR_Egger) == "se"] <- "se_Egger"
      names(MR_Egger)[names(MR_Egger) == "pval"] <- "pval_Egger"
      final_res<- MR_Egger %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_Egger=NA
      final_res$se_Egger=NA
      final_res$pval_Egger=NA
    }
    
    if (exists("IVW")==TRUE){
      names(IVW)[names(IVW) == "b"] <- "b_IVW"
      names(IVW)[names(IVW) == "se"] <- "se_IVW"
      names(IVW)[names(IVW) == "pval"] <- "pval_IVW"
      final_res<- IVW %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_IVW=NA
      final_res$se_IVW=NA
      final_res$pval_IVW=NA
    } 
    
    if (exists("IVW_MRE")==TRUE){
      names(IVW_MRE)[names(IVW_MRE) == "b"] <- "b_IVW_MRE"
      names(IVW_MRE)[names(IVW_MRE) == "se"] <- "se_IVW_MRE"
      names(IVW_MRE)[names(IVW_MRE) == "pval"] <- "pval_IVW_MRE"
      final_res<- IVW_MRE %>% right_join(final_res, by= "id.exposure")
    } else {
      final_res$b_IVW_MRE=NA
      final_res$se_IVW_MRE=NA
      final_res$pval_IVW_MRE=NA
    }
    
    names(final_res)[names(final_res) == "id.exposure"] <- "id"
    
    # final_res<-final_res%>%select(Test, id, trait, b_IVW_MRE, se_IVW_MRE, pval_IVW_MRE, b_Egger, se_Egger, pval_Egger,
    #                               Egger_intercept, pval_intercept, Het_IVW_pval, Het_Egger_pval,
    #                               b_W_Med, se_W_Med, pval_W_Med, b_W_Mod, se_W_Mod, pval_W_Mod,
    #                               b_PRESSO_raw, se_PRESSO_raw, pval_PRESSO_raw, b_PRESSO_corrected, se_PRESSO_corrected, pval_PRESSO_corrected, pval_PRESSO_Global,
    #                               Methods, nsnps, Beta, SE, pval, 
    #                               snp_r2.exposure, snp_r2.outcome, Correct_causal_direction,steiger_pval,
    #                               len_SNP, samplesize.exposure, Total_R_Square, F_stat, F_stat_sim)
    
    final_res<-final_res%>%select(Test, b_IVW_MRE, se_IVW_MRE, pval_IVW_MRE, b_Egger, se_Egger, pval_Egger,
                                  Egger_intercept, pval_intercept, Het_IVW_pval, Het_Egger_pval,
                                  b_W_Med, se_W_Med, pval_W_Med, b_W_Mod, se_W_Mod, pval_W_Mod,
                                  b_PRESSO_raw, se_PRESSO_raw, pval_PRESSO_raw, b_PRESSO_corrected, se_PRESSO_corrected, pval_PRESSO_corrected, pval_PRESSO_Global,
                                  Methods, nsnps, Beta, SE, pval, 
                                  len_SNP, samplesize.exposure, Total_R_Square, F_stat, F_stat_sim)
    
    final_res[is.na(final_res) ] <- "-"
    #####################      @@@@@@@@@@@@@ change file name
    Outputfile=paste(Pathway_out,n, "_WBC.txt", sep="")
    write.table(final_res, file= Outputfile, col.names = FALSE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
  }
}

for (n in Outcome_list) {
  outcomefile=paste(Pathway_SNP,  n ,"_SNP.txt", sep="")
  ####### Change csv file
  exp_dat <- read_exposure_data(
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
  exp_dat$id.exposure = n
  exp_dat$exposure = n
  
  Trait <- read.csv("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/WBC_Trait.txt",header=F, as.is=T,sep = "\t")
  Trait$V4=as.numeric(gsub(",", "", Trait$V4))
  
  len_exp_file=length(Trait$V1)
  for (e in 1:len_exp_file) {
    if (Trait$V2[e]==0 || Trait$V2[e]==10 || Trait$V2[e]==30) {
      outcome_dat <- extract_outcome_data(snps=exp_dat$SNP, outcomes=Trait$V1[e])
      if (length(exp_dat$SNP)>0) {
        My_MR(exp_dat,outcome_dat)
      }
    } else if (Trait$V2[e]==1) {
      exp_dat_file=paste(Pathway_SNP, Trait$V1[e], "_buildGRCh37.tsv", sep="")
      outcome_dat=read_outcome_data(
        filename = exp_dat_file,
        snps = exp_dat$SNP,
        sep = "\t",
        #phenotype_col = "outcome",
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        eaf_col = "effect_allele_frequency",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        #samplesize_col = "all_meta_sample_N",
        pval_col = "p_value")
      outcome_dat$samplesize.outcome = Trait$V4[e]
      outcome_dat$id.outcome = Trait$V1[e]
      outcome_dat$outcome = Trait$V1[e]
      if (length(exp_dat$SNP)>0) {
        My_MR(exp_dat,outcome_dat)
      }
    } else if (Trait$V2[e]==2) {
      exp_dat_file=paste(Pathway_SNP, Trait$V1[e], "_buildGRCh37_out.tsv", sep="")
      outcome_dat=read_outcome_data(
        filename = exp_dat_file,
        snps = exp_dat$SNP,
        sep = "\t",
        #phenotype_col = "outcome",
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        eaf_col = "effect_allele_frequency",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        #samplesize_col = "all_meta_sample_N",
        pval_col = "p_value")
      outcome_dat$samplesize.outcome = Trait$V4[e]
      outcome_dat$id.outcome = Trait$V1[e]
      outcome_dat$outcome = Trait$V1[e]
      if (length(exp_dat$SNP)>0) {
        My_MR(exp_dat,outcome_dat)
      }
    } else if (Trait$V2[e]==3) {
      exp_dat_file=paste(Pathway_SNP, "27863252-", Trait$V1[e], ".tsv", sep="")
      outcome_dat=read_outcome_data(
        filename = exp_dat_file,
        snps = exp_dat$SNP,
        sep = "\t",
        #phenotype_col = "outcome",
        snp_col = "hm_rsid",
        beta_col = "hm_beta",
        se_col = "standard_error",
        eaf_col = "hm_effect_allele_frequency",
        effect_allele_col = "hm_effect_allele",
        other_allele_col = "hm_other_allele",
        #samplesize_col = "all_meta_sample_N",
        pval_col = "p_value")
      outcome_dat$samplesize.outcome = Trait$V4[e]
      outcome_dat$id.outcome = Trait$V1[e]
      outcome_dat$outcome = Trait$V1[e]
      if (length(exp_dat$SNP)>0) {
        My_MR(exp_dat,outcome_dat)
      }
    }
  }
}

