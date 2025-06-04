# Main analysis: Multilevel modeling
library(nlme)
library(ggplot2)
library(reghelper)
library(viridis)
library(boot) 


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


coredf1<- cbind(dataset="DF1",coremm(df.1))
coredf2<- cbind(dataset="DF2",coremm(df.2))
coredf3<- cbind(dataset="DF3",coremm(df.3))
output_coremm<- rbind(
  coredf1,
  coredf2,
  coredf3)
write.csv(output_coremm, paste0("manuscript/results/res_mm_mainanalysis_Table3TableS5_",Sys.Date(),".csv"))

# ===========================================================================
# Sensitivity Analysis: Alternative model specifications (Table S6.2 and S6.3)
# ===========================================================================

# Sensitivity analyses: alternative model specifications
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


outputdf1<- cbind(dataset="DF1",manymm(df.1))
outputdf2<- cbind(dataset="DF2",manymm(df.2))
outputdf3<- cbind(dataset="DF3",manymm(df.3))


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

df.1s1 <- calcSwitching(as.data.frame(dfraw1),inputEmotions1s1)
df.1s2 <- calcSwitching(as.data.frame(dfraw1),inputEmotions1s2)
df.1s3 <- calcSwitching(as.data.frame(dfraw1),inputEmotions1s3)
df.2s1 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s1)
df.2s2 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s2)
df.2s3 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s3)
df.2s4 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2s4)
df.3s1 <- calcSwitching(as.data.frame(dfraw4),inputEmotions3s1)
df.3s2 <- calcSwitching(as.data.frame(dfraw4),inputEmotions3s2)
df.3s3 <- calcSwitching(as.data.frame(dfraw4),inputEmotions3s3)
df.3s4 <- calcSwitching(as.data.frame(dfraw4),inputEmotions3s4)
if(file.exists(filepathEMOTE)){
  df.3s5 <- calcSwitching(as.data.frame(dfraw4),inputEmotions3s5)
  df.3s6 <- calcSwitching(as.data.frame(dfraw4),inputEmotions3s6)
}

LOO_names <- c(paste0("df.1s", 1:3),
               paste0("df.2s", 1:4),
               paste0("df.3s", 1:ifelse(file.exists(filepathEMOTE),6,4))
               )
sensitivityLOO_loop <- do.call(rbind, lapply(LOO_names, function(df_name) {
  df <- get(df_name)      # Get the dataframe
  result <- coremm(df)       # Apply main analysis
  cbind("source" = df_name, # Add the source name as a new column
        result)                 # Return the modified result
}))

write.csv(sensitivityLOO_loop, paste0("manuscript/results/res_mm_sensitivity_LOO_TableS61_",Sys.Date(),".csv"))


# ==================================
# Simple slope analysis
# ==================================

ssmm1 <- mm2(df.1)
ssmm2 <- mm2(df.2)
ssmm3 <- mm2(df.3)

ssgraph2 <- graph_model(ssmm2, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD)

?graph_model
ssdf1 <- simple_slopes(ssmm1)
ssdf2 <- simple_slopes(ssmm2)
ssdf3 <- simple_slopes(ssmm3)

ssdf1.t <- simple_slopes(ssmm1, levels=list(person_CESD=c(0.00,0.28, 'sstest')))
ssdf2.t <- simple_slopes(ssmm2, levels=list(person_CESD=c(0.016,0.183, 'sstest')))
ssdf2.t <- simple_slopes(ssmm2, levels=list(person_CESD=c(0.01,0.19, 'sstest')))
ssdf3.t <- simple_slopes(ssmm3, levels=list(person_CESD=c(0.00,0.13, 'sstest')))
ssdf3.textra <- simple_slopes(ssmm3, levels=list(person_CESD=c(0.115, 'sstest')))
ssdf1.t
ssdf2
ssdf2.t
ssdf3
ssdf3.t

output_ss <- rbind(cbind(source = "DF1", ssdf1),
cbind(source = "DF2", ssdf2),
cbind(source = "DF3", ssdf3))

write.csv(output_ss, paste0("manuscript/results/res_mm_ss_Table4_",Sys.Date(),".csv"))


# <0.00, 0.016, and <0.00 are the corresponding threshold of CES-D in dataset 1, 2, 3
# which participants had to have in order to experience an _increase_ of negative emotion intensity
# This means that according to model predictions, only participants who reported zero CES-D (for all items)
# are expected to experience such an increase
sum(dfperson.2$person_CESD<0.01)/nrow(dfperson.2)
# Which there were none.

# As for those expected to experience a decrease...
(sum(dfperson.1$person_CESD>0.28)+
sum(dfperson.2$person_CESD>0.19)+
sum(dfperson.3$person_CESD>0.13))/(nrow(dfperson.1)+nrow(dfperson.2)+nrow(dfperson.3))


#==============================
# Bootstrapping

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

safeboot.H2 <- function(data, indices) {
  tryCatch(
    {
      # Your original function logic here
      boot.H2(data, indices)
    },
    error = function(e) {
      NA  # Or any default value, or a named list if FUN returns complex objects
    }
  )
}

bootrep <- 1000

set.seed(1999)
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


boot::t
testboot$t0
testboot$t
boot.ci
ci = boot.ci(bootH1df1, type = "perc")
ci
ci = boot.ci(testboot2000, type = "perc")
ci = boot.ci(testboot2000, type = "bca")
bootci_H1df1 <- boot.ci(bootH1df1, type = "bca")
bootci_H1df2 <- boot.ci(bootH1df2, type = "bca")
bootci_H1df3 <- boot.ci(bootH1df3, type = "bca")


#==============================
# Graphs

g1 <- graph_model(mm_CESDcompleteL1D.df1, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD, 
            labels = list("y" = "Intensity of Negative Emotions\n (controlling for previous intensity)",
                          "x" = "Negative Emotion Transition",
                          "lines" = "Level of Baseline\nDepressive Symptoms"),
            errorbars = "none") + theme_bw()
g2 <- graph_model(mm_CESDcompleteL1D.df2, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD, 
                  labels = list("y" = "Intensity of Negative Emotions\n (controlling for previous intensity)",
                                "x" = "Negative Emotion Transition",
                                "lines" = "Level of Baseline\nDepressive Symptoms"),
                  errorbars = "none") + theme_bw()
g3 <- graph_model(mm_CESDcompleteL1D.df3, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD, 
                  labels = list("y" = "Intensity of Negative Emotions\n (controlling for previous intensity)",
                                "x" = "Negative Emotion Transition",
                                "lines" = "Level of Baseline\nDepressive Symptoms"),
                  errorbars = "none") + theme_bw()

g1data <-ggplot_build(g1)
g2data <-ggplot_build(g2)
g3data <-ggplot_build(g3)
N.df1 <- length(unique(mm_CESDcompleteL1D.df1$groups$ppnr))
N.df2 <- length(unique(mm_CESDcompleteL1D.df2$groups$ppnr))
N.df3 <- length(unique(mm_CESDcompleteL1D.df3$groups$ppnr))
gcombine.1.y <- (g1data$data[[1]]$y*N.df1 + g2data$data[[1]]$y*N.df2 + g3data$data[[1]]$y*N.df3) / sum(N.df1,N.df2,N.df3)
# gcombine.1.ymin <- (g1data$data[[1]]$ymin*N.df1 + g2data$data[[1]]$ymin*N.df2 + g3data$data[[1]]$ymin*N.df3) / sum(N.df1,N.df2,N.df3)
gcombine.1.ymax <- (g1data$data[[1]]$ymax*N.df1 + g2data$data[[1]]$ymax*N.df2 + g3data$data[[1]]$ymax*N.df3) / sum(N.df1,N.df2,N.df3)

gcombine.2.y <- (g1data$data[[2]]$y*N.df1 + g2data$data[[2]]$y*N.df2 + g3data$data[[2]]$y*N.df3) / sum(N.df1,N.df2,N.df3)
# gcombine.2.ymin <- (g1data$data[[2]]$ymin*N.df1 + g2data$data[[2]]$ymin*N.df2 + g3data$data[[2]]$ymin*N.df3) / sum(N.df1,N.df2,N.df3)
gcombine.2.ymax <- (g1data$data[[2]]$ymax*N.df1 + g2data$data[[2]]$ymax*N.df2 + g3data$data[[2]]$ymax*N.df3) / sum(N.df1,N.df2,N.df3)

gcombine_raw <- graph_model(mm_CESDcompleteL1D.df2, y=moment_NA, x=moment_NA_bray.bal.succw, lines=person_CESD, 
                  labels = list("y" = "Negative Emotion Intensity",
                                "x" = "Negative Emotion Transition",
                                "lines" = "Level of Baseline\nDepressive Symptoms"),
                  errorbars = "none", ymin = min(gcombine.1.ymax), ymax = max(gcombine.1.ymax)) + theme_bw()
gcombinedata <- ggplot_build(gcombine_raw)
gcombinedata$data[[1]]$y <- gcombine.1.y
gcombinedata$data[[1]]$ymax <- gcombine.1.ymax
gcombinedata$data[[2]]$y <- gcombine.2.y
gcombinedata$data[[2]]$ymax <- gcombine.2.ymax
gcombinedata$layout$panel_params[[1]]$y.range

gnew <- ggplot_gtable(gcombinedata)
grid::grid.newpage()
grid::grid.draw(gnew) 


varlist_graph <- c("ppnr","moment_NA_bray.bal.suc","moment_NA","moment_NAcwL1D","timecw","person_CESD")
df.graph <- rbind.data.frame(
  cbind.data.frame("study" = "DF1",df.1[,varlist_graph]),
  cbind.data.frame("study" = "DF2",df.2[,varlist_graph]),
  cbind.data.frame("study" = "DF3",df.3[,varlist_graph]))
modelgraph <- lme(fixed=moment_NA ~ (moment_NA_bray.bal.suc)*person_CESD +timecw+
                    + moment_NAcwL1D,
                  data=df.graph,
                  random=~1+ (moment_NA_bray.bal.suc)*person_CESD+ moment_NAcwL1D | ppnr, correlation = corAR1(),
                  control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
modelgraph <- lme(fixed=moment_NA ~ (moment_NA_bray.bal.suc)*person_CESD +timecw,
                  data=df.graph,
                  random=~1+ (moment_NA_bray.bal.suc)*person_CESD | ppnr, correlation = corAR1(),
                  control =list(msMaxIter = 1000, msMaxEval = 1000, opt = "optim"),na.action = na.omit)
summary(modelgraph)

model1 <- mm_simple.df1
model1 <- mm_CESDcompleteL1D.df1

#obtaining predicted scores for individuals
df.graph$rownum <- seq.int(nrow(df.graph))
df.graph.pred <- data.frame(rownum = (rownames(modelgraph$residuals)),
           pred = predict(modelgraph))
df.graph.proto <- data.frame(rownum = (rownames(modelgraph$residuals)),
                         proto = predict(modelgraph, level=0:1))
df.graph <- merge(x = df.graph, y = df.graph.pred, by = "rownum", all = TRUE)
df.graph <- merge(x = df.graph, y = df.graph.proto, by = "rownum", all = TRUE)

df.graph$bin_CESD <- findInterval(df.graph$person_CESD, (0:10)/10)
df.graph <- df.graph[order(df.graph$person_CESD, decreasing = TRUE), ]
df.graph <- df.graph[order(df.graph$person_CESD, ), ]

ordered_ids <- df.graph %>%
  group_by(ppnr) %>%
  summarize(mean_value = mean(person_CESD)) %>%
  arrange(mean_value) %>%    # ascending: first drawn, last on top
  pull(ppnr)

#plotting predicted trajectories
ggplot(data = df.graph, aes(x = moment_NA_bray.bal.suc, y = pred, group = ppnr, color=person_CESD)) +
  # geom_point() +
  # geom_line(alpha=.2) +
  geom_line(aes(x = moment_NA_bray.bal.suc, y = proto.predict.ppnr),size=1)+
  # lapply(ordered_ids, function(i) {
  #   geom_line(data = df.graph[df.graph$ppnr == i, ],
  #             aes(x = moment_NA_bray.bal.suc, y = proto), size = 1)
  # }) +
    labs(color = "Baseline\nDepressive\nSymptoms")+
  xlab("Negative Emotion Transition (t-1 to t)") + 
  ylab("Negative Emotion Intensity (t)\nControlling for Intensity (t-1)") + #ylim(0,60) +
  scale_x_continuous(breaks=seq(0,1,by=0.2)) + 
  scale_color_gradient2(midpoint=0.5, low="#F8E621", mid="#1F938B", #F9E3D1, #C51751, #35193E, midpoint=mean(df.graph$person_CESD),
                       high="#450558", space ="Lab" , limits = c(0, 1))+
  theme_bw()
?scale_color_gradient2
#obtaining predicted scores for individuals
df.1$rownum <- seq.int(nrow(df.1))
df.1.pred <- data.frame(rownum = (rownames(model1$residuals)),
                        pred.model1 = predict(model1))
df.1.proto <- data.frame(rownum = (rownames(model1$residuals)),
                         proto.model1 = predict(model1, level=0))
df.1 <- merge(x = df.1, y = df.1.pred, by = "rownum", all = TRUE)
df.1 <- merge(x = df.1, y = df.1.proto, by = "rownum", all = TRUE)

df.1$bin_CESD <- findInterval(df.1$person_CESD, (0:10)/10)
#plotting predicted trajectories
ggplot(data = df.1, aes(x = moment_NA_bray.bal.succw, y = pred.model1, group = ppnr, color=person_CESD)) +
  # geom_point() +
  # geom_line(alpha=.2) +
  geom_line(aes(x = moment_NA_bray.bal.succw, y = proto.model1),size=1) +
  labs(color = "Baseline\nDepressive\nSymptoms")+
  xlab("Negative Emotion Transition") + 
  ylab("Negative Emotion Intensity") + ylim(0,100) +
  scale_x_continuous(breaks=seq(-0.5,1,by=0.5)) + scale_color_gradient2(midpoint=mean(df.1$person_CESD), low="blue", mid="grey",
                                                                        high="red", space ="Lab" )+
  theme_bw()

