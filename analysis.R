# Main analysis: Multilevel modeling
# Only run this after running "prepare_data.R"!

library(nlme)
library(ggplot2)
library(reghelper)
library(viridis)
library(boot) 
library(ggpubr)


# ===========================================================================
# Main Analysis: Model 1 and Model 2 in their original specifications
#     Corresponds to Table 3 and Table S5.
# ===========================================================================


# Testing the main hypotheses (H1 & H2) with Model 1 and Model 2

mm1 <- function (df){
  lme(fixed=moment_NA ~ moment_NA_bray.bal.succw+moment_NA_bray.gra.succw+moment_NAcwL1D+timecw+
        moment_NA_bray.bal.succb+moment_NA_bray.gra.succb,
      data=df,
      random=~1+ moment_NA_bray.bal.succw+moment_NA_bray.gra.succw+moment_NAcwL1D | ppnr, correlation = corAR1(),
      control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
}
mm2 <- function (df){
  lme(fixed=moment_NA ~ (moment_NA_bray.bal.succw+moment_NA_bray.gra.succw)*person_CESD +timecw+
        moment_NA_bray.bal.succb+moment_NA_bray.gra.succb + moment_NAcwL1D,
      data=df,
      random=~1+ (moment_NA_bray.bal.succw+moment_NA_bray.gra.succw)*person_CESD+ moment_NAcwL1D | ppnr, correlation = corAR1(),
      control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
}

# A wrapper function to prepare multilevel model output

preparemmresult <- function (m){
  cbind(model = deparse(substitute(m)),
        FEest = summary(m)$tTable[,1],
        SE = summary(m)$tTable[,2],
        DF = summary(m)$tTable[,3],
        pvalue = summary(m)$tTable[,5],
        residual = m$sigma^2,
        phi = coef(m$modelStruct$corStruct, unconstrained = FALSE), # needs coef to make this work 
        nobs = nobs(m),
        N = m$dims$ngrps[1],
        u95CI = summary(m)$tTable[,1]+summary(m)$tTable[,2]*1.96,
        l95CI = summary(m)$tTable[,1]-summary(m)$tTable[,2]*1.96,
        AIC = summary(m)$AIC,
        BIC = summary(m)$BIC,
        LL = summary(m)$logLik,
        RMSE = sqrt(mean(m$residuals^2))
  )
  
}

# A wrapper function to output main analysis results

coremm <- function(dfmm){
  
  # Hypothesis 1
  mm_H1 <- mm1(dfmm)
  mm_H2 <- mm2(dfmm)
  
  rbind(
    preparemmresult(mm_H1),
    preparemmresult(mm_H2)
  )
}

# Run main analyses for the 3 datasets
coredf1<- cbind(dataset="DF1",coremm(df.1))
coredf2<- cbind(dataset="DF2",coremm(df.2))
coredf3<- cbind(dataset="DF3",coremm(df.3))

# Bind them and output
output_coremm<- rbind(
  coredf1,
  coredf2,
  coredf3)
write.csv(output_coremm, paste0("manuscript/results/res_mm_mainanalysis_Table3TableS5_",Sys.Date(),".csv"))

# ===========================================================================
# Sensitivity Analysis: Alternative model specifications (Table S6.2 and S6.3)
# ===========================================================================

# Wrapper for sensitivity analyses: alternative model specifications
manymm <- function(dfmm, skipb = FALSE){
  
  # Hypothesis 1 - non-centered dissimilarity
  mm_H1r <- lme(fixed=moment_NA ~ moment_NA_bray.bal.suc+moment_NA_bray.gra.suc+moment_NAcwL1D+timecw,
                data=dfmm,
                random=~1+ moment_NA_bray.bal.suc+moment_NA_bray.gra.suc+moment_NAcwL1D | ppnr, correlation = corAR1(),
                control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
  # Hypothesis 2 - non-centered dissimilarity
  mm_H2r <- lme(fixed=moment_NA ~ (moment_NA_bray.bal.suc+moment_NA_bray.gra.suc)*person_CESD +timecw+
                  + moment_NAcwL1D,
                data=dfmm,
                random=~1+ (moment_NA_bray.bal.suc+moment_NA_bray.gra.suc)*person_CESD+ moment_NAcwL1D | ppnr, correlation = corAR1(),
                control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
  # Hypothesis 1 - without nestedness
  mm_H1s <- lme(fixed=moment_NA ~ moment_NA_bray.bal.succw+moment_NAcwL1D+timecw+
                  moment_NA_bray.bal.succb,
                data=dfmm,
                random=~1+ moment_NA_bray.bal.succw+moment_NAcwL1D | ppnr, correlation = corAR1(),
                control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
  # Hypothesis 2 - without nestedness
  mm_H2s <- lme(fixed=moment_NA ~ (moment_NA_bray.bal.succw)*person_CESD +timecw+
                  moment_NA_bray.bal.succb + moment_NAcwL1D,
                data=dfmm,
                random=~1+ (moment_NA_bray.bal.succw)*person_CESD+ moment_NAcwL1D | ppnr, correlation = corAR1(),
                control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
  rbind(
    preparemmresult(mm_H1s),
    preparemmresult(mm_H1r),
    preparemmresult(mm_H2s),
    preparemmresult(mm_H2r)
  )
}

# Run sensitivity analyses for the 3 datasets
outputdf1<- cbind(dataset="DF1",manymm(df.1))
outputdf2<- cbind(dataset="DF2",manymm(df.2))
outputdf3<- cbind(dataset="DF3",manymm(df.3))

# Bind them and output
output_manymm<- rbind(
  outputdf1,
  outputdf2,
  outputdf3)
write.csv(output_manymm, paste0("manuscript/results/res_mm_sensitivity_TableS62andS63_",Sys.Date(),".csv"))



#============================
# Sensitivity Analysis - Leave one out (LOO)
#============================

# Leave-one-out specifications on which emotions to include to derive Bray-Curtis dissimilarity

inputEmotions1s1 <- c("M_nieder","M_bekuem")
inputEmotions1s2 <- c("M_nerv","M_bekuem")
inputEmotions1s3 <- c("M_nerv","M_nieder")
inputEmotions2s1 <- c("SAD","ANX","DEPRE")
inputEmotions2s2 <- c("ANGRY", "ANX","DEPRE")
inputEmotions2s3 <- c("ANGRY", "SAD","DEPRE")
inputEmotions2s4 <- c("ANGRY", "SAD","ANX")

if(file.exists(filepathEMOTE)){
  inputEmotions3s1 <- c("droevig","angstig","depressief","stress","lonely")
  inputEmotions3s2 <- c("kwaad","angstig","depressief","stress","lonely")
  inputEmotions3s3 <- c("kwaad","droevig","depressief","stress","lonely")
  inputEmotions3s4 <- c("kwaad","droevig","angstig","stress","lonely")
  inputEmotions3s5 <- c("kwaad","droevig","angstig","depressief","lonely")
  inputEmotions3s6 <- c("kwaad","droevig","angstig","depressief","stress")
}else{
  inputEmotions3s1 <- c("droevig","angstig","depressief")
  inputEmotions3s2 <- c("kwaad","angstig","depressief")
  inputEmotions3s3 <- c("kwaad","droevig","depressief")
  inputEmotions3s4 <- c("kwaad","droevig","angstig")
}

# Re-calculating emotion transitions with the new sets of emotion items
df.1s1 <- calcSwitching(as.data.frame(dfraw1),inputEmotions1s1)
df.1s2 <- calcSwitching(as.data.frame(dfraw1),inputEmotions1s2)
df.1s3 <- calcSwitching(as.data.frame(dfraw1),inputEmotions1s3)
df.2s1 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s1)
df.2s2 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s2)
df.2s3 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s3)
df.2s4 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s4)
df.3s1 <- calcSwitching(as.data.frame(dfraw3),inputEmotions3s1)
df.3s2 <- calcSwitching(as.data.frame(dfraw3),inputEmotions3s2)
df.3s3 <- calcSwitching(as.data.frame(dfraw3),inputEmotions3s3)
df.3s4 <- calcSwitching(as.data.frame(dfraw3),inputEmotions3s4)
if(file.exists(filepathEMOTE)){
  df.3s5 <- calcSwitching(as.data.frame(dfraw3),inputEmotions3s5)
  df.3s6 <- calcSwitching(as.data.frame(dfraw3),inputEmotions3s6)
}


# Define the list of leave one out analysis to run
LOO_names <- c(paste0("df.1s", 1:3),
               paste0("df.2s", 1:4),
               paste0("df.3s", 1:ifelse(file.exists(filepathEMOTE),6,4))
               )

# A loop function that runs all those LOO analyses
# Warning! It might take an hour to run the below.
sensitivityLOO_loop <- do.call(rbind, lapply(LOO_names, function(df_name) {
  df <- get(df_name)      # Get the dataframe
  result <- coremm(df)       # Apply main analysis
  cbind("source" = df_name, # Add the source name as a new column
        result)                 # Return the modified result
}))

# Output
write.csv(sensitivityLOO_loop, paste0("manuscript/results/res_mm_sensitivity_LOO_TableS61_",Sys.Date(),".csv"))

# ==================================
# Simple slope analysis
# ==================================

# Obtain multilevel model results for Model 2
ssmm1 <- mm2(df.1)
ssmm2 <- mm2(df.2)
ssmm3 <- mm2(df.3)

# Simple slope analyses for the 3 datasets
ssdf1 <- simple_slopes(ssmm1)
ssdf2 <- simple_slopes(ssmm2)
ssdf3 <- simple_slopes(ssmm3)

# JARS requires 95%CI  
ssdf1$l95ci <- ssdf1$`Test Estimate`-1.96*ssdf1$`Std. Error`
ssdf1$u95ci <- ssdf1$`Test Estimate`+1.96*ssdf1$`Std. Error`
ssdf2$l95ci <- ssdf2$`Test Estimate`-1.96*ssdf2$`Std. Error`
ssdf2$u95ci <- ssdf2$`Test Estimate`+1.96*ssdf2$`Std. Error`
ssdf3$l95ci <- ssdf3$`Test Estimate`-1.96*ssdf3$`Std. Error`
ssdf3$u95ci <- ssdf3$`Test Estimate`+1.96*ssdf3$`Std. Error`

# Finding out the thresholds of CES-D (according to model predictions)
# where participants start to experience a significant (p <.05) transition-intensity reduction effect
ssdf1.t <- simple_slopes(ssmm1, levels=list(person_CESD=c(0.00,0.28, 'sstest')))
ssdf2.t <- simple_slopes(ssmm2, levels=list(person_CESD=c(0.01,0.19, 'sstest')))
ssdf3.t <- simple_slopes(ssmm3, levels=list(person_CESD=c(0.00,0.115, 'sstest')))

# bind results and output
output_ss <- rbind(cbind(source = "DF1", ssdf1),
cbind(source = "DF2", ssdf2),
cbind(source = "DF3", ssdf3))

write.csv(output_ss, paste0("manuscript/results/res_mm_ss_Table4_",Sys.Date(),".csv"))


# <0.00, 0.016, and <0.00 are the corresponding threshold of CES-D in dataset 1, 2, 3
# which participants had to have in order to experience an _increase_ of negative emotion intensity
# This means that according to model predictions, only participants who reported zero CES-D (for all items)
# are expected to experience such an increase
sum(dfperson.2$person_CESD<0.016)/nrow(dfperson.2)
# Which there were none.

# As for those expected to experience a decrease...
# Bind results and output
output_proportion <- data.frame(
          overall = (sum(dfperson.1$person_CESD>0.28)+
                      sum(dfperson.2$person_CESD>0.19)+
                      sum(dfperson.3$person_CESD>0.13))/(nrow(dfperson.1)+nrow(dfperson.2)+nrow(dfperson.3)),
        dataset1 = sum(dfperson.1$person_CESD>0.28)/nrow(dfperson.1),
        dataset2 = sum(dfperson.2$person_CESD>0.19)/nrow(dfperson.2),
        dataset3 = sum(dfperson.3$person_CESD>0.13)/nrow(dfperson.3))

write.csv(output_proportion, paste0("manuscript/results/res_mm_ss_proportion.csv"))


#==============================
# Bootstrapping
# Warning! Adjust the reps before you run it or it can take like days before R stops.
#==============================

# Warning! It takes quite long to run 1000 repetitions (around 2 days on taktsun's laptop)
bootrep <- 1000

boot.H1 <- function(data, indices) {
  val <- data[indices,] # selecting sample with boot
  fit <-  mm1(val)
  return(fit$coefficients$fixed["moment_NA_bray.bal.succw"])
}

boot.H2 <- function(data, indices) {
  val <- data[indices,] # selecting sample with boot
  fit <- mm2(val)
  return(fit$coefficients$fixed["moment_NA_bray.bal.succw:person_CESD"])
}

# we set this because in some resampling occasions in bootstrapping, there were nlme evaluation erros
safeboot.H2 <- function(data, indices) {
  tryCatch(
    {
      boot.H2(data, indices)
    },
    error = function(e) {
      NA  
    }
  )
}


# set seed to get reproducible results
set.seed(1999)

# we save the boot results in case R crashes during the long bootstrapping process.

bootH1df1 <- boot(data=df.1, statistic=boot.H1,
                  R=bootrep, parallel = "snow")
saveRDS(bootH1df1, paste0("temp/boot",bootrep,"_df1_H1.RData"))
bootH1df2 <- boot(data=df.2, statistic=boot.H1,
                  R=bootrep, parallel = "snow")
saveRDS(bootH1df2, paste0("temp/boot",bootrep,"_df2_H1.RData"))
bootH1df3 <- boot(data=df.3, statistic=boot.H1,
                  R=bootrep, parallel = "snow")
saveRDS(bootH1df3, paste0("temp/boot",bootrep,"_df3_H1.RData"))

bootH2df1 <- boot(data=df.1, statistic=safeboot.H2,
                  R=bootrep, parallel = "snow")
saveRDS(bootH2df1, paste0("temp/boot",bootrep,"_df1_H2.RData"))
bootH2df2 <- boot(data=df.2, statistic=safeboot.H2,
                  R=bootrep, parallel = "snow")
saveRDS(bootH2df2, paste0("temp/boot",bootrep,"_df2_H2.RData"))
bootH2df3 <- boot(data=df.3, statistic=safeboot.H2,
                  R=bootrep, parallel = "snow")
saveRDS(bootH2df3, paste0("temp/boot",bootrep,"_df3_H2.RData"))


output_boot <- rbind(cbind(dataset = "DF1", h="H1",boot.ci(bootH1df1, type = "perc")$percent),
      cbind(dataset = "DF2", h="H1",boot.ci(bootH1df2, type = "perc")$percent),
      cbind(dataset = "DF3", h="H1",boot.ci(bootH1df3, type = "perc")$percent),
      cbind(dataset = "DF1", h="H2",boot.ci(bootH2df1, type = "perc")$percent),
      cbind(dataset = "DF2", h="H2",boot.ci(bootH2df2, type = "perc")$percent),
      cbind(dataset = "DF3", h="H2",boot.ci(bootH2df3, type = "perc")$percent)
)

citest <- boot.ci(bootH1df1, type = "all")
write.csv(output_boot,"manuscript/results/res_bootstrap_TableS5.csv")

#==============================
# Simple Slope Graphs (Figure 3)
#==============================


outssg1 <- graph_model(ssmm1, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD, 
            labels = list("title" = "Dataset 1",
                          "x" = "",
                          "y" = "Overall Negative Emotions Intensity (t)\n Controlling for the Previous Intensity (t-1)",
                          "lines" = "Level of Baseline\nDepressive Symptoms"),
            errorbars = "none", ymin=0, ymax=35) + theme_bw()  + scale_color_brewer(palette = "Set2")

outssg2 <- graph_model(ssmm2, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD, 
            labels = list("title" = "Dataset 2",
                          "y" = " ",
                          "x" = "Negative Emotion Transition (t-1 to t)",
                          "lines" = "Level of Baseline\nDepressive Symptoms"),
            errorbars = "none", ymin=0, ymax=35) + theme_bw() + scale_color_brewer(palette = "Set2")

outssg3 <- graph_model(ssmm3, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD, 
            labels = list("title" = "Dataset 3",
                          "y" = " ",
                          "x" = " ",
                          "lines" = "Level of Baseline\nDepressive Symptoms"),
            errorbars = "none", ymin=0, ymax=35) + theme_bw() + scale_color_brewer(palette = "Set2")
ggarrange(outssg1, outssg2, outssg3, ncol = 3, nrow = 1, common.legend = TRUE, legend="right")

