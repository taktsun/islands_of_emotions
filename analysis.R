# Main analysis: Multilevel modeling
# Only run this after running "prepare_data.R"!

# INSTALL PACKAGES IF NEEDED
# install.packages("devtools")
# devtools::install_github("seanchrismurphy/emodiff")


library(nlme)
library(ggplot2)
library(reghelper)
library(viridis)
library(boot)
library(ggpubr)
library(viridis) 
library(colorspace)
library(emodiff)

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
# Sensitivity Analysis: Alternative model specifications (Table S7.1 and S7.2)
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
write.csv(output_manymm, paste0("manuscript/results/res_mm_sensitivity_TableS71andS72_",Sys.Date(),".csv"))



#============================
# Sensitivity Analysis: Leave one out (LOO) (Table S6)
#============================

# Leave-one-out specifications on which emotions to include to derive Bray-Curtis dissimilarity

inputEmotions1s1 <- c("M_nieder","M_bekuem") # exclude nervous
inputEmotions1s2 <- c("M_nerv","M_bekuem") # exclude
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
write.csv(sensitivityLOO_loop, paste0("manuscript/results/res_mm_sensitivity_LOO_TableS6_",Sys.Date(),".csv"))


#==================================
# Sensitivity Analysis: LOO Graph
#==================================

output_coreLOO <- rbind(output_coremm, sensitivityLOO_loop)
dfoutcoreLOO <- as.data.frame(output_coreLOO)
dfoutcoreLOO$parameter <- rownames(output_coreLOO)


# Rename dataset names in a consistent manner for easy sorting later
dfoutcoreLOO <- dfoutcoreLOO %>%
  mutate(dataset = recode(dataset,
                          "DF1" = "df.1s0",
                          "DF2" = "df.2s0",
                          "DF3" = "df.3s0"))

# Create two strings that map dataset values and the text we want to display
(strLOOmap1 <- unique(dfoutcoreLOO$dataset))
strLOOmap2 <- c("Dataset 1 (complete)",
                "Dataset 2 (complete)",
                "Dataset 3 (complete)",
                "  Excl. Nervous",
                "  Excl. Downhearted",
                "  Excl. Distressed",
                "  Excl. Angry",
                "  Excl. Sad",
                "  Excl. Anxious",
                "  Excl. Depressed",
                "  Excl. Angry",
                "  Excl. Sad",
                "  Excl. Anxious",
                "  Excl. Depressed")       

# Handle the case when other researchers do not have access to the full EMOTE database
if (file.exists(filepathEMOTE)) strLOOmap2 <- c(strLOOmap2,
                                                "  Excl. Stressed",
                                                "  Excl. Lonely")
# Map these labels to the dataframe
dfoutcoreLOO$label <- strLOOmap2[ match(dfoutcoreLOO$dataset, strLOOmap1) ]


# Prepare extra empty rows so that there will be spaces between datasets
extra_rows <- data.frame(dataset = c("df.0s0","df.1s9", "df.2s9"),
                         label = c("","",""),
                         # replicate NA for all other columns in df_plot
                         lapply(df_plot[ , !(names(df_plot) %in% c("dataset","label"))], 
                                function(x) NA))


# Plot for Hypothesis 1

df_plot <- dfoutcoreLOO %>%
  filter(parameter == "moment_NA_bray.bal.succw",
         model     == "mm_H1") %>%
  # Ensure numeric CI/estimate columns
  mutate(across(c(FEest, l95CI, u95CI), as.numeric)) %>%
  # In case l/u are swapped in any row
  mutate(l95 = pmin(l95CI, u95CI),
         u95 = pmax(l95CI, u95CI)) %>%
  # Lock a single, consistent y mapping for BOTH geoms
  mutate(dataset = factor(dataset, levels = unique(dataset)))
graphlabelLOO <- "Hypothesis 1: Negative Emotion Transitions \nPredicting Changes in Overall Emotion Intensity"

# Bind them to the existing dataframe
df_plot <- rbind(df_plot, extra_rows)
df_plot <- df_plot %>%
  mutate(
    # extract 4th and 6th character
    ds        = as.numeric(substr(dataset, 4, 4)), # dataset
    condition = as.numeric(substr(dataset, 6, 6)) # LOO condition
  ) %>%
  # order dataset factor by ds first, then condition
  arrange(ds, condition) %>%
  mutate(dataset = factor(dataset, levels = rev(unique(dataset))))
df_plot$color <- df_plot$ds*100 + df_plot$condition*10+20
df_plot$line_size <- ifelse(df_plot$condition == 0, 1.5, 1)
df_plot$point_size <- ifelse(df_plot$condition == 0, 3, 1.8)



# Create plot.
# The three rows of warnings are expected because we purposely created empty rows to separate the 3 datasets in the graph
ggplot(df_plot, aes(y = dataset, x = FEest,
                    colour = color)) +
  geom_segment(aes(x = l95, xend = u95, yend = dataset),
               linewidth = df_plot$line_size, lineend = "round") +
  geom_point(size = df_plot$point_size) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_y_discrete(labels = setNames(df_plot$label, df_plot$dataset)) +
  labs(x = "Fixed Effect Estimates (Dots) and 95% CI (Lines)", y = "") +
  theme_classic() + scale_colour_viridis_b()+
  theme(
    legend.position = "none",        # remove legend
    axis.ticks.y = element_blank(),  # remove tick marks
    axis.line.y  = element_blank(),   # remove y-axis line
    axis.text.y = element_text(hjust=0)
  )+
  labs(title = graphlabelLOO)

# Create plot for hypothesis 2
df_plot <- dfoutcoreLOO %>%
  filter(parameter == "moment_NA_bray.bal.succw:person_CESD",
         model     == "mm_H2") %>%
  # Ensure numeric CI/estimate columns
  mutate(across(c(FEest, l95CI, u95CI), as.numeric)) %>%
  # In case l/u are swapped in any row
  mutate(l95 = pmin(l95CI, u95CI),
         u95 = pmax(l95CI, u95CI)) %>%
  # Lock a single, consistent y mapping for BOTH geoms
  mutate(dataset = factor(dataset, levels = unique(dataset)))
graphlabelLOO <- "Hypothesis 2: Depressive Symptoms Moderate The Within-Person\nAssociation Between Emotion Transition and Intensity Change"

# Bind them to the existing dataframe
df_plot <- rbind(df_plot, extra_rows)
df_plot <- df_plot %>%
  mutate(
    # extract 4th and 6th character
    ds        = as.numeric(substr(dataset, 4, 4)),
    condition = as.numeric(substr(dataset, 6, 6))
  ) %>%
  # order dataset factor by ds first, then condition
  arrange(ds, condition) %>%
  mutate(dataset = factor(dataset, levels = rev(unique(dataset))))
df_plot$color <- df_plot$ds*100 + df_plot$condition*10+20
df_plot$line_size <- ifelse(df_plot$condition == 0, 1.5, 1)
df_plot$point_size <- ifelse(df_plot$condition == 0, 3, 1.8)



# Create plot.
# The three rows of warnings are expected because we purposely created empty rows to separate the 3 datasets in the graph
ggplot(df_plot, aes(y = dataset, x = FEest,
                    colour = color)) +
  geom_segment(aes(x = l95, xend = u95, yend = dataset),
               linewidth = df_plot$line_size, lineend = "round") +
  geom_point(size = df_plot$point_size) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_y_discrete(labels = setNames(df_plot$label, df_plot$dataset)) +
  labs(x = "Fixed Effect Estimates (Dots) and 95% CI (Lines)", y = "") +
  theme_classic() + scale_colour_viridis_b()+
  theme(
    legend.position = "none",        # remove legend
    axis.ticks.y = element_blank(),  # remove tick marks
    axis.line.y  = element_blank(),   # remove y-axis line
    axis.text.y = element_text(hjust=0)
  )+
  labs(title = graphlabelLOO)


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


# ==========================================
# Marginal effects graph
# ============================================

# Define level of depressive symptoms: 
# 0 (min), 1 (max),
# M -2, -1, 0, +1, +2 SD
# Levels at which the moderator starts to be significant

ssvalue1 <- c(0,
              mean(dfperson.1$person_CESD),
              mean(dfperson.1$person_CESD) - 2*sd(dfperson.1$person_CESD),
              mean(dfperson.1$person_CESD) - 1*sd(dfperson.1$person_CESD),
              mean(dfperson.1$person_CESD) + 1*sd(dfperson.1$person_CESD),
              mean(dfperson.1$person_CESD) + 2*sd(dfperson.1$person_CESD),
              mean(dfperson.1$person_CESD) + 3*sd(dfperson.1$person_CESD),
              0.27,
              0.28,
              1
)

# Sort the values in ascending order
ssvalue1 <- sort(ssvalue1)
# Remove impossible values (negative values or >1) 
ssvalue1 <- ssvalue1[ssvalue1 >= 0 & ssvalue1 <= 1]

ssvalue2 <- c(0,
              mean(dfperson.2$person_CESD),
              mean(dfperson.2$person_CESD) - 2*sd(dfperson.2$person_CESD),
              mean(dfperson.2$person_CESD) - 1*sd(dfperson.2$person_CESD),
              mean(dfperson.2$person_CESD) + 1*sd(dfperson.2$person_CESD),
              mean(dfperson.2$person_CESD) + 2*sd(dfperson.2$person_CESD),
              mean(dfperson.2$person_CESD) + 3*sd(dfperson.2$person_CESD),
              0.01,
              0.18,
              0.19,
              1
)
ssvalue2 <- sort(ssvalue2)
ssvalue2 <- ssvalue2[ssvalue2 >= 0 & ssvalue2 <= 1]

ssvalue3 <- c(0,
              mean(dfperson.3$person_CESD),
              mean(dfperson.3$person_CESD) - 2*sd(dfperson.3$person_CESD),
              mean(dfperson.3$person_CESD) - 1*sd(dfperson.3$person_CESD),
              mean(dfperson.3$person_CESD) + 1*sd(dfperson.3$person_CESD),
              mean(dfperson.3$person_CESD) + 2*sd(dfperson.3$person_CESD),
              mean(dfperson.3$person_CESD) + 3*sd(dfperson.3$person_CESD),
              0.11,
              0.12,
              1
)
ssvalue3 <- sort(ssvalue3)
ssvalue3 <- ssvalue3[ssvalue3 >= 0 & ssvalue3 <= 1]

# Calculate the marginal effects & 95% CI at different levels of depressive symptoms
ssdf1.tr <- simple_slopes(ssmm1, levels=list(person_CESD=ssvalue1), confint = TRUE)
ssdf2.tr <- simple_slopes(ssmm2, levels=list(person_CESD=ssvalue2), confint = TRUE)
ssdf3.tr <- simple_slopes(ssmm3, levels=list(person_CESD=ssvalue3), confint = TRUE)

# Extract relevant columns (variables) for further graph plotting
ros1 <- ssdf1.tr[,2:6]
ros2 <- ssdf2.tr[,2:6]
ros3 <- ssdf3.tr[,2:6]

# Wrapper for producing the graph
graph_marginal <- function(strxlab = "Depressive Symptoms (CES-D: Scaled 0-1)",
                           strylab = "Marginal Effect of Negative Emotion Transition on Intensity Change"){

  # Create a temp df that stores the frequencies of different values of depressive symptoms  
  df_rug <- df.ros |>
    count(person_CESD, name = "freq")
  
  # Change the column names to that used by McCabe et al. (2018)'s interActive shinyapp
  colnames(ros) <- c("hyp.Z","pe.X","SE","LL","UL")

  # The below scripts are adapted from McCabe et al. (2018)'s interActive shinyapp
  
  # Determine whether marginal effects are significant or not
  ros$significance[(ros$LL*ros$UL) > 0] <- "sig"
  ros$significance[(ros$LL*ros$UL) <= 0] <- "not sig"
  # Determine the sign of effect (positive or negative)
  ros$sign[ros$pe.X<0]<-"neg"
  ros$sign[ros$pe.X>=0]<-"pos"
  
  # Specify the moderator variable name
  mod = "person_CESD"
  # Ensure things are in numeric format so that the script will work
  ros$hyp.Z        <- as.numeric(ros$hyp.Z)
  
  #Creating the RoS plot
  rosplot <- ggplot() +
    geom_ribbon(data = ros, aes(x = hyp.Z, ymin = LL, ymax = UL),
                fill=brewer.pal(3,"Greys")[3], alpha = .25) +
    geom_line(data = ros, aes(hyp.Z, pe.X), color="Black", size = 1.25) +
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1),limits = c(0,1)) +
    ylim(min(-45),max(12)) +
    geom_hline(yintercept = 0) +
    xlab(strxlab) +
    ylab(paste(strylab)) +
    geom_rug(
      data = df_rug,
      aes(x = person_CESD, size = freq),
      inherit.aes = FALSE,
      sides = "b"
    ) +
    scale_size(range = c(0.2, 2), guide = "none")+
    theme_bw() +
    theme(text=element_text(family="Helvetica",size=12, color="black"))
  
  # Draw vertical dashed lines to indicate the threshold levels at which marginal effects start to be significant
  if("not sig"%in%ros[which(ros$sign=="neg"),"significance"]==TRUE &
     "sig"%in%ros[which(ros$sign=="neg"),"significance"]==TRUE ){  
    
    sigline.negneg<-min(ros[which(ros$significance=="sig" & ros$sign=="neg"),"hyp.Z"])
    rosplot <- rosplot + 
      geom_vline(xintercept = sigline.negneg, linetype="longdash")
  }
  
  if("not sig"%in%ros[which(ros$sign=="pos"),"significance"]==TRUE &
     "sig"%in%ros[which(ros$sign=="pos"),"significance"]==TRUE){  
    
    sigline.negpos<-max(ros[which(ros$significance=="sig" & ros$sign=="pos"),"hyp.Z"])
    rosplot <- rosplot + 
      geom_vline(xintercept = sigline.negpos, linetype="dotdash")
  }
  # Draw the plot
  rosplot
}

ros <- ros1
df.ros<-dfperson.1
graph_marginal(strxlab = "") # graph for study1 doesn't need x axis
ros <- ros2
df.ros<-dfperson.2
graph_marginal(strylab = "") # graph for study 2 doesn't need a y-axis
ros <- ros3
df.ros<-dfperson.3
graph_marginal(strxlab = "",strylab = "") # graph for study 3 doesn't need either



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



#===========================================
# Sensitivity Analysis: Emotion differentiation as covariates
#===========================================


df.1ED <- emodiff::calculate_ed(dat = df.1[df.1$b_completeNA,], emotions = inputEmotions1, ppnr, allow_neg_icc = TRUE)
df.1ED$m_EDL1 <- lagvar(m_ED, id=ppnr, obs=triggerid, data=df.1ED)
df.1ED$person_ED <- calc.mean(m_ED, ppnr, data=df.1ED,expand=TRUE)
df.1ED$grand_m_ED <- mean(calc.mean(m_ED, ppnr, data=df.1ED),na.rm=TRUE)
df.1ED$m_EDcb <- df.1ED$person_ED - df.1ED$grand_m_ED


df.2ED <- emodiff::calculate_ed(dat = df.2[df.2$b_completeNA,], emotions = inputEmotions2, ppnr, allow_neg_icc = TRUE)
df.2ED$m_EDL1 <- lagvar(m_ED, id=ppnr, obs=triggerid, data=df.2ED)
df.2ED$person_ED <- calc.mean(m_ED, ppnr, data=df.2ED,expand=TRUE)
df.2ED$grand_m_ED <- mean(calc.mean(m_ED, ppnr, data=df.2ED),na.rm=TRUE)
df.2ED$m_EDcb <- df.2ED$person_ED - df.2ED$grand_m_ED


df.3ED <- emodiff::calculate_ed(dat = df.3[df.3$b_completeNA,], emotions = inputEmotions3, ppnr, allow_neg_icc = TRUE)
df.3ED$m_EDL1 <- lagvar(m_ED, id=ppnr, obs=triggerid, data=df.3ED)
df.3ED$person_ED <- calc.mean(m_ED, ppnr, data=df.3ED,expand=TRUE)
df.3ED$grand_m_ED <- mean(calc.mean(m_ED, ppnr, data=df.3ED),na.rm=TRUE)
df.3ED$m_EDcb <- df.3ED$person_ED - df.3ED$grand_m_ED


mm1mcED <- function (df){
  lme(fixed=moment_NA ~ moment_NA_bray.bal.succw+moment_NA_bray.gra.succw+moment_NAcwL1D+m_EDL1+timecw+
        moment_NA_bray.bal.succb+moment_NA_bray.gra.succb+m_EDcb,
      data=df,
      random=~1+ moment_NA_bray.bal.succw+moment_NA_bray.gra.succw+moment_NAcwL1D+m_EDL1 | ppnr, correlation = corAR1(),
      control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
}
mm2mcED <- function (df){
  lme(fixed=moment_NA ~ (moment_NA_bray.bal.succw+moment_NA_bray.gra.succw)*person_CESD +timecw+
        moment_NA_bray.bal.succb+moment_NA_bray.gra.succb + moment_NAcwL1D+m_EDL1+m_EDcb,
      data=df,
      random=~1+ (moment_NA_bray.bal.succw+moment_NA_bray.gra.succw)*person_CESD+ moment_NAcwL1D+m_EDL1 | ppnr, correlation = corAR1(),
      control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
}



manyEDmm <- function(dfmm){
  rbind(
    preparemmresult(mm1mcED(dfmm)),
    preparemmresult(mm2mcED(dfmm))
  )
}

# Run main analyses for the 3 datasets
manyEDdf1<- cbind(dataset="DF1",manyEDmm(df.1ED))
manyEDdf2<- cbind(dataset="DF2",manyEDmm(df.2ED))
manyEDdf3<- cbind(dataset="DF3",manyEDmm(df.3ED))

# Bind them and output
output_manyEDmm<- rbind(
  manyEDdf1,
  manyEDdf2,
  manyEDdf3)
write.csv(output_manyEDmm, paste0("manuscript/results/res_mm_sensitivity_TableS73_",Sys.Date(),".csv"))


