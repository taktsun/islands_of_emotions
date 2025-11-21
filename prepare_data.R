# In this file, write the R-code necessary to load your original data file
# (e.g., an SPSS, Excel, or SAS-file), and convert it to a data.frame. Then,
# use the function open_data(your_data_frame) or closed_data(your_data_frame)
# to store the data.

# you may need to install the esmpack package with the below 2 lines
# install.packages("remotes")
# remotes::install_github("wviechtb/esmpack")

library(worcs)
library(RCurl)
library(haven) # for reading sav files
library(lubridate) #for function day
library(dplyr)
library(esmpack) # for the lagvar function
library(betapart)

filepathEMOTE <- 'temp/data_downloads_C7SC6HWU8R_2025-04-01Leuven_3-wave_longitudinal_study.csv'

# ===========================================================================
# List the negative emotion items we are interested in
# ===========================================================================

inputEmotions1 <- c("M_nerv","M_nieder","M_bekuem")
inputEmotions2 <- c("ANGRY",
                    "SAD",
                    "ANX",
                    "DEPRE")
if (file.exists(filepathEMOTE)){
  inputEmotions3 <- c("kwaad","droevig","angstig","depressief","stress","lonely")
}else{
  inputEmotions3 <- c("kwaad","droevig","angstig","depressief") # see read raw data explanatory note below
}

# ===========================================================================
# Read raw data and the dissimilarity calculation script
# (read locally if not connected to the Internet)
# ===========================================================================

# Note! The public version of Dataset 3 (raw3) on OSF contains only 4 negative emotions.
# However, the full dataset contains 6 negative emotions.
# We have requested the 6 negative emotion items from EMOTE under data request code "C7SC6HWU8R"
# with which you can request the same dataset from EMOTE (https://emotedatabase.com/)

# A wrapper to check whether you have an internet connection 
has_internet <- function(url = "https://osf.io") {
  tryCatch({
    getURL(url, timeout = 5)
    TRUE
  }, error = function(e) FALSE)
}

bInternet <- has_internet("github.com")
if (bInternet){
  source("https://github.com/taktsun/dissimilarity-for-ESM-data/raw/main/BrayCurtisDissimilarity_Calculate.R")
}else{
  source("temp/BrayCurtisDissimilarity_Calculate.R")
}

bInternet <- has_internet("osf.io")
if (bInternet){
  # be careful about time-out errors!
  # If needed, download these datasets to your local drive to load them...
  raw1 <- read_sav('https://osf.io/download/w8y33/')
  raw2 <- read_sav('https://osf.io/download/gm52c/')
  raw3 <- read_sav('https://osf.io/download/uvqjh/')
}else{
  # locally stored copies of data
  raw1 <- read_sav('temp/BData1.sav')
  raw2 <- read_sav('temp/BData2.sav')
  raw3 <- read_sav('temp/BData3.sav')
}
if (file.exists(filepathEMOTE)){
  raw3 <- read.csv(filepathEMOTE)
  raw3 <- raw3[raw3$WAVE_ES==1,]
  # Named vector: names are old, values are new
  raw3newnames <- c(ANG_ES = "kwaad", SAD_ES = "droevig", FEAR_ES = "angstig", DEP_ES = "depressief", STR_ES = "stress", LONE_ES = "lonely")
  # Rename matching columns
  names(raw3)[names(raw3) %in% names(raw3newnames)] <- raw3newnames[names(raw3)[names(raw3) %in% names(raw3newnames)]]

}



# ===========================================================================
# Functions for data pre-processing
# ===========================================================================

# detect whether there are zero variance across multiple variables across ESM observations (pre-registered exclusion criterion)
# from: https://github.com/taktsun/ED_ERV/blob/master/func_preprocessing.R
zerovariance <- function(dftemp,variables_to_check){
  # Group by PPID and check if the specified variables remain the same in all rows
  dftemp <- dftemp[complete.cases(dftemp[, variables_to_check]), ]
  same_values_df <- dftemp %>%
    group_by(ppnr) %>%
    filter(all(across(all_of(variables_to_check), ~. == first(.))))
  # Extract the list of PPID values
  unique_ppid_list <- unique(same_values_df$ppnr)
  # Print the list of PPID values
  unique_ppid_list
}


prepRawData <- function(rawdata, sID, sBeep, sTS, sDateFormat, cEmotion, cCESD){
  d <- data.frame(rawdata)
  colnames(d)[which(names(d) == sID)] <- "ppnr"
  colnames(d)[which(names(d) == sBeep)] <- "triggerid"
  colnames(d)[which(names(d) == sTS)] <- "timestamp"
  d$timestamp <- as_datetime(d$timestamp)
  dateformat<-sDateFormat
  d <- d %>%
    group_by(ppnr) %>%
    mutate(
      earliest_time = min(timestamp, na.rm = TRUE),
      timesince = as.numeric(difftime(timestamp, earliest_time, units = "secs"))
    ) %>%
    select(-earliest_time) %>%  # optional: remove helper column
    ungroup()
  d$day <-yday(parse_date_time(d$timestamp,orders=dateformat))
  d <- d %>%
    group_by(ppnr) %>%
    mutate(day = day - min(day) + 1) %>%
    ungroup()
  inputEmotions <- cEmotion
  inputCESD <- cCESD
  inputNeeded <- append(c("ppnr","triggerid","day","timesince"),inputEmotions)
  d <- d[,inputNeeded]
  d$person_CESD <- rowMeans(rawdata[,inputCESD],na.rm = TRUE)
  
  d
}

prepRawDataEMOTE <- function(rawdata, sID, sBeep, cEmotion, cCESD, sTimesince, sDay){
  d <- data.frame(rawdata)
  colnames(d)[which(names(d) == sID)] <- "ppnr"
  colnames(d)[which(names(d) == sBeep)] <- "triggerid"
  colnames(d)[which(names(d) == sTimesince)] <- "timesince"
  colnames(d)[which(names(d) == sDay)] <- "day"
  d <- d %>%
    group_by(ppnr) %>%
    mutate(day = day - min(day) + 1) %>%
    ungroup()
  inputEmotions <- cEmotion
  inputCESD <- cCESD
  inputNeeded <- append(c("ppnr","triggerid","day","timesince"),inputEmotions)
  d <- d[,inputNeeded]
  d$person_CESD <- rowMeans(rawdata[,inputCESD],na.rm = TRUE)
  
  d
}


# ===========================================================================
# Data pre-processing: reverse item, harmonize naming of variables across 3 datasets
# ===========================================================================


# Three datasets have different ranges of CESD. Harmonization needed later.
range(raw1[,paste0("T1_CESD",c(1:20))])
range(raw2[,paste0("cesd",c(1:20))],na.rm=TRUE)
if(file.exists(filepathEMOTE)){
       range(raw3[,paste0("CESD_",c(1:20),"_BL")], na.rm = TRUE)
  }else{
       range(raw3[,paste0("cesd",c(1:20),"_w1")])
  }


# Reverse item handling
raw1$T1_CESD4 <- 4 - raw1$T1_CESD4
raw1$T1_CESD8 <- 4 - raw1$T1_CESD8
raw1$T1_CESD12 <- 4 - raw1$T1_CESD12
raw1$T1_CESD16 <- 4 - raw1$T1_CESD16
raw2$cesd4 <- 3 - raw2$cesd4
raw2$cesd8 <- 3 - raw2$cesd8
raw2$cesd12 <- 3 - raw2$cesd12
raw2$cesd16 <- 3 - raw2$cesd16
if(file.exists(filepathEMOTE)){
  raw3$CESD_4_BL <- 3 - raw3$CESD_4_BL
  raw3$CESD_8_BL <- 3 - raw3$CESD_8_BL
  raw3$CESD_12_BL <- 3 - raw3$CESD_12_BL
  raw3$CESD_16_BL <- 3 - raw3$CESD_16_BL
}else{

  raw3$cesd4_w1 <- 3 - raw3$cesd4_w1
  raw3$cesd8_w1 <- 3 - raw3$cesd8_w1
  raw3$cesd12_w1 <- 3 - raw3$cesd12_w1
  raw3$cesd16_w1 <- 3 - raw3$cesd16_w1
  }

# Extract just the variables we need...

dfraw1 <- prepRawData(rawdata = raw1,
                      sID = "ID_anonym",
                      sBeep = "a_ftl_0",
                      sTS = "A_time",
                      sDateFormat = "ymd HMS",
                      cEmotion = inputEmotions1,
                      cCESD = paste0("T1_CESD",seq(1:20)))

dfraw2 <- prepRawData(rawdata = raw2,
                      sID = "PpID",
                      sBeep = "BeepNo",
                      sTS = "date_time",
                      sDateFormat = "ymd HMS",
                      cEmotion = inputEmotions2,
                      cCESD = paste0("cesd",seq(1:20)))
if(file.exists(filepathEMOTE)){
  dfraw3 <- raw3
  # convert date into integer with lubridate::dmy
  dfraw3$intDate <- as.integer(dmy(dfraw3$Date_Local))
  # minus the integer with the minimum value within each participant
  # and then +1, so that each participant starts from DAY=1
  dfraw3 <- dfraw3 %>%
    group_by(UUID) %>%
    mutate(DAY = intDate - min(intDate) + 1)
  # convert timestamp to integer (in seconds of a 24-hour day)
  # and plus DAY*86400 so that we can calculate beeps
  dfraw3$intSec <- as.integer(hms(dfraw3$Time_Local))+dfraw3$DAY*86400
  dfraw3 <- dfraw3 %>%
    group_by(UUID) %>%
    mutate(BEEPTIME = intSec - min(intSec) + 1)
  dfraw3 <- dfraw3 %>%
    group_by(UUID) %>%
    mutate(BEEP = order(BEEPTIME))
  dfraw3$BEEPTIME <- dfraw3$BEEPTIME-1
  dfraw3 <- prepRawDataEMOTE(rawdata = dfraw3,
                        sID = "UUID",
                        sBeep = "BEEP",
                        cEmotion = c(inputEmotions3),
                        cCESD = paste0("CESD_",seq(1:20),"_BL"),
                        sTimesince = "BEEPTIME",
                        sDay = "DAY")
  dfraw3$ppnr <- as.numeric(factor(dfraw3$ppnr, labels = 1:length(unique(dfraw3$ppnr))))
}else{
  dfraw3 <- prepRawData(rawdata = raw3,
                        sID = "PpID",
                        sBeep = "BeepNr",
                        sTS = "date_time",
                        sDateFormat = "ymd HMS",
                        cEmotion = c(inputEmotions3),
                        cCESD = paste0("cesd",seq(1:20),"_w1"))
}

# ===========================================================================
# Data pre-processing: Harmonizing the range items
#     so that all emotion items range from 0 to 100
#     and all person CESD scores range from 0 to 1
# ===========================================================================

# As noted earlier, harmonize CESD to a range from 0 to 1
dfraw1$person_CESD <- dfraw1$person_CESD/4
dfraw2$person_CESD <- dfraw2$person_CESD/3
dfraw3$person_CESD <- dfraw3$person_CESD/3

# Three datasets have different scales for ESM items
range(dfraw1[,inputEmotions1],na.rm=TRUE)
range(dfraw2[,inputEmotions2],na.rm=TRUE)
range(dfraw3[,inputEmotions3],na.rm=TRUE)


# Therefore, harmonize the rating of negative emotion items
dfraw1[,inputEmotions1]<-(dfraw1[,inputEmotions1]*100/6) 
dfraw2[,inputEmotions2]<-(dfraw2[,inputEmotions2]-1)/99*100 


# They now have the same range
range(dfraw1[,inputEmotions1],na.rm=TRUE)
range(dfraw2[,inputEmotions2],na.rm=TRUE)
range(dfraw3[,inputEmotions3],na.rm=TRUE)

# ===========================================================================
# Excluding participants based on potentially problematic responses
# ===========================================================================

# Did any participants have zero variance across all negative emotion items?
# nobody showed zero variance, so no participants need to be excluded
zerovariance(dfraw1,inputEmotions1)
zerovariance(dfraw2,inputEmotions2)
zerovariance(dfraw3,inputEmotions3)


# Did any participants give the same rating for all CESD items?
# Note: CESD has reverse items. So same ratings for all CESD items meantindicate problematic responses.

# Dataset1: no participants had probleatic CESD
dfraw1$CESDequal <- apply(raw1[,paste0("T1_CESD",seq(1:20))], 1, function(row) all(row == row[1]))
unique(dfraw1$ppnr[dfraw1$CESDequal])
# Dataset2: ppid 82 has problematic CESD
dfraw2$CESDequal <- apply(raw2[,paste0("cesd",seq(1:20))], 1, function(row) all(row == row[1]))
unique(dfraw2$ppnr[dfraw2$CESDequal])
# Dataset3: ppid 170 has problematic CESD
if(file.exists(filepathEMOTE)){
  dfraw3$CESDequal <- apply(raw3[,paste0("CESD_",seq(1:20),"_BL")], 1, function(row) all(row == row[1]))
  }else{
    dfraw3$CESDequal <- apply(raw3[,paste0("cesd",seq(1:20),"_w1")], 1, function(row) all(row == row[1]))
}
unique(dfraw3$ppnr[dfraw3$CESDequal])


# Based on the above, we 1 participant from Dataset 2 and 1 participant from Dataset 3
dfraw2 <- dfraw2[dfraw2$ppnr != unique(dfraw2$ppnr[dfraw2$CESDequal]),]
dfraw3 <- dfraw3[dfraw3$ppnr != unique(dfraw3$ppnr[dfraw3$CESDequal]),]

# ===========================================================================
# Calculate Bray-Curtis dissimilarity
# ===========================================================================

calcSwitching <- function(df, inputEmotions){
  #identify the first beep within each person
  df$firstbeep <- df$triggerid==ave(df$triggerid, df$ppnr, FUN = min)
  # boolean: does the observation has complete NA ratings (i.e., no missingness)?
  df$b_completeNA <- complete.cases(df[,inputEmotions])
  if(min(df$triggerid,na.rm = TRUE)==0){
    df$triggerid <- df$triggerid +1
  }
  # Create boolean of which rows suitable for successive comparisons. Three conditions:
  # 1. if the observation is NOT the first beep of a person
  # 2. if the moment of interest has complete ER strategy ratings (i.e., no missingness)
  # 3. if the moment before the moment of interest has complete ER strategy ratings
  
  
  df$b_completeNAL1 <- lagvar(b_completeNA, id=ppnr, obs=triggerid, data=df)
  df$b_completeNAL1[is.na(df$b_completeNAL1)] <- FALSE
  
  
  df$b_successive <- !df$firstbeep & df$b_completeNA & df$b_completeNAL1
  
  
  
  # calculate successive dissimilarity
  
  df <- calcBrayCurtisESM(df, inputEmotions, "ppnr", "triggerid")
  df$moment_NA_bray.all.suc <- df$BrayCurtisFull.suc   # if one of two rows are all-zero, return 1
  df$moment_NA_bray.bal.suc <- df$BrayCurtisRepl.suc   # if one of two rows are all-zero, return NA
  df$moment_NA_bray.gra.suc <- df$BrayCurtisNest.suc   # if one of two rows are all-zero, return NA
  
  
  df$timecw <- calc.mcent(triggerid, ppnr, data=df)
  
  df$moment_NA <- rowMeans(df[,c(inputEmotions)], na.rm=TRUE)
  df$moment_NAcw <- calc.mcent(moment_NA, ppnr, data=df)
  df$moment_NAL1D <- lagvar(moment_NA,ppnr,obs=triggerid,day=day, data = df)
  df$moment_NAcwL1D <- lagvar(moment_NAcw,ppnr,obs=triggerid,day=day, data = df)
  df$person_NA <- calc.mean(moment_NA, ppnr, data=df,expand=TRUE)
  
  
  # calculate moment-level variability indices centered within a person
  df$moment_NA_bray.all.succw <-calc.mcent(moment_NA_bray.all.suc, ppnr, data=df)
  df$moment_NA_bray.bal.succw <-calc.mcent(moment_NA_bray.bal.suc, ppnr, data=df)
  df$moment_NA_bray.gra.succw <-calc.mcent(moment_NA_bray.gra.suc, ppnr, data=df)
  
  # calculate person-level variability indices,
  # i.e., value is the same across all time points within a person
  df$person_NA_bray.all.suc <- calc.mean(moment_NA_bray.all.suc, ppnr, data=df,expand=TRUE)
  df$person_NA_bray.bal.suc <- calc.mean(moment_NA_bray.bal.suc, ppnr, data=df,expand=TRUE)
  df$person_NA_bray.gra.suc <- calc.mean(moment_NA_bray.gra.suc, ppnr, data=df,expand=TRUE)
  
  # calculate grand means
  df$grand_moment_NA <- mean(calc.mean(moment_NA, ppnr, data=df),na.rm=TRUE)
  df$grand_moment_NA_bray.all.suc <- mean(calc.mean(moment_NA_bray.all.suc, ppnr, data=df),na.rm=TRUE)
  df$grand_moment_NA_bray.bal.suc <- mean(calc.mean(moment_NA_bray.bal.suc, ppnr, data=df),na.rm=TRUE)
  df$grand_moment_NA_bray.gra.suc <- mean(calc.mean(moment_NA_bray.gra.suc, ppnr, data=df),na.rm=TRUE)
  df$grand_person_CESD <- mean(calc.mean(person_CESD, ppnr, data=df),na.rm=TRUE)
  
  
  # calculate between component
  df$moment_NAcb <- df$person_NA - df$grand_moment_NA
  df$moment_NA_bray.all.succb <- df$person_NA_bray.all.suc - df$grand_moment_NA_bray.all.suc
  df$moment_NA_bray.bal.succb <- df$person_NA_bray.bal.suc - df$grand_moment_NA_bray.bal.suc
  df$moment_NA_bray.gra.succb <- df$person_NA_bray.gra.suc - df$grand_moment_NA_bray.gra.suc
  
  
  df
  
}

df.1 <- calcSwitching(as.data.frame(dfraw1),inputEmotions1)
df.2 <- calcSwitching(as.data.frame(dfraw2),inputEmotions2)
df.3 <- calcSwitching(as.data.frame(dfraw3),inputEmotions3)

