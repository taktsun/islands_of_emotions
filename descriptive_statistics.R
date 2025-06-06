library(dplyr)
library(misty) # for calculating ICC in descriptive statistics
library(psych) # multilevel correlation
library(multilevelTools) # for omega reliability

#============================
# Descriptive statistics: Age, Gender, and how many people are above CES-D cutoffs?
#============================

# Descriptives extracted from previous studies
# because these data aren't available in the public dataset

# Dataset 1: Extracted from Blanke et al. (2020):  http://dx.doi.org/10.1037/emo0000566
agemean_df1 <- 25.55
agesd_df1 <- 2.74
N_df1 <- 70
female_df1 <- 35/N_df1
agemin_df1 <- 20
agemax_df1 <- 30

# Dataset 2: Extracted from Blanke et al. (2020):  http://dx.doi.org/10.1037/emo0000566
agemean_df2 <- 19.06
agesd_df2 <- 1.28
N_df2 <- 95
female_df2 <- 59/N_df2
agemin_df2 <- 18
agemax_df2 <- 24

# Dataset 3: Extracted from Lo et al. (2025): https://doi.org/10.1007/s42761-025-00301-4 (supplemental materials)
agemean_df3 <- 18.32
agesd_df3 <- 0.96
N_df3 <- 202
female_df3 <- 0.55
agemin_df3 <- 17
agemax_df3 <- 24


ageinfo <- data.frame( min = min(agemin_df1,agemin_df2,agemin_df3),
                       
                      max = max(agemax_df1,agemax_df2,agemax_df3),
                      mean = (agemean_df1*N_df1 + agemean_df2*N_df2 + agemean_df3*N_df3)/sum(N_df1,N_df2,N_df3),
                      sd = sqrt((agesd_df1^2*(N_df1-1) + agesd_df2^2*(N_df2-1) + agesd_df3^2*(N_df3-1))/sum(N_df1,N_df2,N_df3-3)),
                      female = (female_df1*N_df1 + female_df2*N_df2 + female_df3*N_df3)/sum(N_df1,N_df2,N_df3))

write.csv(ageinfo, "manuscript/results/methods_participants_age.csv")

# Distribution of CES-D

dfperson.1 <- distinct(df.1, ppnr, .keep_all = TRUE)
dfperson.2 <- distinct(df.2, ppnr, .keep_all = TRUE)
dfperson.3 <- distinct(df.3, ppnr, .keep_all = TRUE)

outputReliabilityCESD <- data.frame(omegadf1 = omega(raw1[!duplicated(raw1$ID_anonym),paste0("T1_CESD",c(1:20))])$omega.tot,
           omegadf2 = omega(raw2[!duplicated(raw2$PpID),paste0("cesd",c(1:20))])$omega.tot,
           omegadf3 = if(file.exists(filepathEMOTE)){
  omega(raw3[!duplicated(raw3$UUID),paste0("CESD_",c(1:20),"_BL")])$omega.tot
}else{
  omega(raw3[!duplicated(raw3$PpID),paste0("cesd",c(1:20),"_w1")])$omega.tot
})

write.csv(outputReliabilityCESD,"manuscript/results/methods_measures_omegaReliability_CESD.csv")
# Cutoff recommended: Tradionally: 16; Meta-analysis: 20
sum(dfperson.1$person_CESD > (16/60))/nrow(dfperson.1)
sum(dfperson.2$person_CESD > (16/60))/nrow(dfperson.2)
sum(dfperson.3$person_CESD > (16/60))/nrow(dfperson.3)

resDesStat_aboveCESDcutoff <- data.frame(dataset1 = sum(dfperson.1$person_CESD > (20/60))/nrow(dfperson.1),
dataset2 = sum(dfperson.2$person_CESD > (20/60))/nrow(dfperson.2),
dataset3 = sum(dfperson.3$person_CESD > (20/60))/nrow(dfperson.3))
resDesStat_aboveCESDcutoff
write.csv(resDesStat_aboveCESDcutoff,"manuscript/results/desstat_aboveCESDcutoff.csv")


#============================
# Descriptive statistics: ESM compliance rate
#============================


resDesStat_complance <- (sum(df.1$b_completeNA)/nrow(df.1)+
    sum(df.2$b_completeNA)/nrow(df.2)+
    sum(df.3$b_completeNA)/nrow(df.3))/3
write.csv(resDesStat_complance, "manuscript/results/desstat_compliance.csv")

#============================
# Descriptive statistics: ESM compliance rate
#============================

# Wrapper functions

desstat <- function(x,group, nameoverride = ""){
  c(
    # variable name
    var = ifelse(nameoverride == "", deparse(substitute(x)),nameoverride),
    # n obs
    nobs = sum(!is.na(x)),
    # grand mean
    mean = round(mean(x, na.rm = TRUE),2),
    # within person SD
    wSD = round(mean(unlist(aggregate(x, by=list(group), FUN=sd, na.rm = TRUE)[2]),na.rm=TRUE),2),
    # between person SD
    bSD = round(sd(unlist(aggregate(x, by=list(group), FUN=mean, na.rm = TRUE)[2]),na.rm=TRUE),2),
    # ICC: how much of the variance is between-person
    ICC = round(multilevel.icc(x, cluster = group),2),
    min =  round(mean(unlist(aggregate(x, by=list(group), FUN=min, na.rm = TRUE)[2]),na.rm=TRUE),2),
    max =  round(mean(unlist(aggregate(x, by=list(group), FUN=max, na.rm = TRUE)[2]),na.rm=TRUE),2)
  )
}


summarydesstat <- function(df){
  rbind(
    desstat(df$moment_NA,df$ppnr),
    desstat(df$moment_NA_bray.all.suc,df$ppnr),
    desstat(df$moment_NA_bray.bal.suc,df$ppnr),
    desstat(df$moment_NA_bray.gra.suc,df$ppnr),
    c("df$person_CESD",nrow(df[df$b_firstbeep,]),round(mean(df[df$b_firstbeep,"person_CESD"]),2),
      0, round(sd(df[df$b_firstbeep,"person_CESD"]),2), 0,
      round(min(df[df$b_firstbeep,"person_CESD"]),2),
      round(max(df[df$b_firstbeep,"person_CESD"]),2))
  )
}

dynamicdesstate <- function(input, df){
  desstat(df[,input],df$ppnr, paste0(substitute(df),"_",input))
}



desstat1 <- summarydesstat(df.1)
desstat2 <- summarydesstat(df.2)
desstat3 <- summarydesstat(df.3)
desstat_all <- rbind((cbind(dataset = 1, desstat1)),
                     (cbind(dataset = 2, desstat2)),
                     (cbind(dataset = 3, desstat3)))
desstat_item <- rbind(
  do.call(rbind, lapply(inputEmotions1, dynamicdesstate, df.1)),
  do.call(rbind, lapply(inputEmotions2, dynamicdesstate, df.2)),
  do.call(rbind, lapply(inputEmotions3, dynamicdesstate, df.3))
  
)
write.csv(desstat_item,"manuscript/results/desstat_TableS4_part1.csv")


#============================
# Multilevel correlation
#============================

# 3 wrapper functions for handling matrix triangles
combineTriangles <- function(matrix1, matrix2) {
  # Create a new matrix with the same dimensions as the input matrices
  new_matrix <- matrix(0, nrow = nrow(matrix1), ncol = ncol(matrix1))
  # Fill in the lower triangle with values from matrix1
  new_matrix[lower.tri(new_matrix)] <- matrix1[lower.tri(matrix1)]
  # Fill in the upper triangle with values from matrix2
  new_matrix[upper.tri(new_matrix)] <- matrix2[upper.tri(matrix2)]
  return(new_matrix)
}
matrixFromTriangles <- function(mlength, diagonal_values, lv, uv) {
  mat <- matrix(0, nrow = mlength, ncol = mlength)
  mat[lower.tri(mat, diag = FALSE)] <- uv
  mat <- t(mat)
  mat[lower.tri(mat, diag = FALSE)] <- lv
  return(mat)
}
cortablesummary <- function(df, inputIndices = inputIndices, labelIndices = labelIndices){
  objStatsby <- statsBy((df[, c("ppnr",inputIndices)]),
                        "ppnr", cors=TRUE)
  tmp.r <- combineTriangles(objStatsby$rwg, objStatsby$rbg) #lower, upper
  tmp.low.ci <- matrixFromTriangles(length(inputIndices),0,
                                    objStatsby$ci.wg$r.ci$lower,objStatsby$ci.bg$r.ci$lower)
  tmp.up.ci <- matrixFromTriangles(length(inputIndices),0,
                                   objStatsby$ci.wg$r.ci$upper,objStatsby$ci.bg$r.ci$upper)
  mat <- matrix(paste0(sub("0.", ".",format(round(tmp.r,2),nsmall=2)), " [",
                       sub("0.", ".",format(round(tmp.low.ci,2),nsmall=2)),",",
                       sub("0.", ".",format(round(tmp.up.ci,2),nsmall=2)),"]"), nrow = nrow(tmp.r),
                dimnames = list(paste0(c(1:length(labelIndices)),". ",labelIndices),
                                c(1:length(labelIndices)))
  )
  
  diag(mat) <- ""
  mat
}
cortablesummaryblank <- function(labelIndices = labelIndices){
  mat <- matrix(0, nrow = length(labelIndices), ncol = (length(labelIndices)),
                dimnames = list(paste0(c(1:length(labelIndices)),". ",labelIndices),
                                c(1:length(labelIndices))))
  mat <- cbind(matrix(0, nrow = length(labelIndices), ncol = 6),
               mat)
  colnames(mat) <- c("n",
                     "M",
                     "SDw",
                     "SDb",
                     "Min",
                     "Max",
                     colnames(mat)[-c(1:6)])
  mat
}

inputCommon <- c("moment_NA","moment_NA_bray.all.suc","moment_NA_bray.bal.suc","moment_NA_bray.gra.suc" ,
                 "person_CESD")

# There are warnings because person level CESD has no within-person SD 
cortablesummary(df.1,inputCommon,inputCommon)
cortablesummary(df.2,inputCommon,inputCommon)
cortablesummary(df.3,inputCommon,inputCommon)

desstat_all <- rbind((cbind(dataset = 1, desstat1,cortablesummary(df.1,inputCommon,inputCommon))),
                     (cbind(dataset = 2, desstat2,cortablesummary(df.2,inputCommon,inputCommon))),
                     (cbind(dataset = 3, desstat3,cortablesummary(df.3,inputCommon,inputCommon))))

# multilevel correlation for individual emotion items.
# Some young adults show zero variance in some emotion items; 4 occasions in total.
desstat_itemdf1 <- cbind(do.call(rbind, lapply(inputEmotions1, dynamicdesstate, df.1)),
                         cortablesummary(df.1,inputEmotions1,inputEmotions1))
desstat_itemdf2 <- cbind(do.call(rbind, lapply(inputEmotions2, dynamicdesstate, df.2)),
                         cortablesummary(df.2,inputEmotions2,inputEmotions2))
desstat_itemdf3 <- cbind(do.call(rbind, lapply(inputEmotions3, dynamicdesstate, df.3)),
                         cortablesummary(df.3,inputEmotions3,inputEmotions3))
write.csv(desstat_all,"manuscript/results/desstat_Table2.csv")
write.csv(desstat_itemdf1,"manuscript/results/desstat_TableS4_df1.csv")
write.csv(desstat_itemdf2,"manuscript/results/desstat_TableS4_df2.csv")
write.csv(desstat_itemdf3,"manuscript/results/desstat_TableS4_df3.csv")

#============================
# Omega Reliability
#============================

# There are certain warnings about zero variance in single ESM measures, 
# which is expected, because we only excluded participants who show zero variance
# in the same GROUP of variables (e.g., in all negative emotion items), but not in distinct items
omega_NA_SM_1 <- omegaSEM(items = inputEmotions1,  data = df.1,  id = "ppnr",  savemodel = FALSE)
omega_NA_SM_2 <- omegaSEM(items = inputEmotions2,  data = df.2,  id = "ppnr",  savemodel = FALSE)
omega_NA_SM_3 <- omegaSEM(items = inputEmotions3,  data = df.3,  id = "ppnr",  savemodel = FALSE)

outputReliability <- rbind(omega_NA_SM_1$Results,
                           omega_NA_SM_2$Results,
                           omega_NA_SM_3$Results)
write.csv(outputReliability, "manuscript/results/methods_measures_omegaReliability_ESM.csv", row.names=FALSE)

