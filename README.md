
# Readme <a href='https://osf.io/zcvbs/'><img src='worcs_icon.png' align="right" height="139" /></a>

This is the github of my manuscript "Negative Emotion Transitions May Have Immediate Benefits in Decreasing Negative Emotions in Daily Life". You can reproduce the analysis results in the manuscript following this readme.

## Data availability

Our study analyzes 3 open datasets. Two of them are hosted on OSF and one hosted on EMOTE:

Dataset 1: https://osf.io/download/w8y33/
Dataset 2: https://osf.io/download/gm52c/
Dataset 3 (Full items): http://emotedatabase.com/requestid/C7SC6HWU8R (Request code: C7SC6HWU8R; parent dataset: https://emotedatabase.com/datasets/22/)

Requesting Dataset 3 from EMOTE takes a few days. In the mean time, it is possible to analyze an alternative version hosted on OSF:
Dataset 3 (Incomplete items): https://osf.io/download/uvqjh/
Note, however, this OSF version only has 4 negative emotion items, instead of 6 in the full dataset from EMOTE. So, the results will be slightly different.

## Where do I start?

You can load this project in RStudio by opening the file called 'EmotionTransition.Rproj'.
If you want to fully reproduce our results, please download Dataset 3 from EMOTE (Request code: C7SC6HWU8R)
and put it under the "temp" folder.

## Project structure

### Folders

 - The manuscript folder has a README.md file that points readers to the latest available pre-print.  
 - The manuscript/results folder is needed as it is the directory where the output files are saved. 
 - The temp folder is empty on Github. Please put Dataset 3 downloaded from EMOTE (data_downloads_C7SC6HWU8R_2025-04-01Leuven_3-wave_longitudinal_study.csv) to this folder.

### R scripts
File                      | Description                      | Usage         
------------------------- | -------------------------------- | --------------
prepare_data.R            | Script to process raw data       | Run to reproduce results
descriptive_statistics.R            | Script to produce descriptive statistics       | Run to reproduce results
analysis.R            | Script to produce output from main and supplemental analyses       | Run to reproduce results

### Project files
File                      | Description                      | Usage         
------------------------- | -------------------------------- | --------------
README.md                 | Description of project           | Read only
EmotionTransition.Rproj   | Project file                     | Loads project 
LICENSE                   | User permissions                 | Read only     
.worcs                    | WORCS metadata YAML              | Read only     
manuscript/README.md | Contains a URL to the latest pre-print            | Read only
renv.lock                 | Reproducible R environment       | Read only     


# Reproducibility

Reproduce the results by these steps.

 1. Install RStudio and R
 2. Install WORCS dependencies
		
		install.packages("worcs", dependencies = TRUE)
		tinytex::install_tinytex()
		renv::consent(provided = TRUE)
		
 3. [Clone](https://resources.github.com/github-and-rstudio/#:~:text=Clone%20the%20repository%20with%20RStudio&text=On%20GitHub%2C%20navigate%20to%20the,RStudio%20on%20your%20local%20environment.) this repo (https://github.com/taktsun/islands_of_emotions) to your RStudio
 4. Restore the package dependencies
	

	    renv::restore()
	    
	    
 5. Download the raw dataset from EMOTE (http://emotedatabase.com/requestid/C7SC6HWU8R) to folder temp. 
 6. Run 3 R scripts in the below sequence to reproduce the results.
 
	- prepare_data.R 
	- descriptive_statistics.R 
	- analysis.R 
All the scripts has accommodated the (un)availability of Dataset 3 from EMOTE. If it doesn't detect the dataset (data_downloads_C7SC6HWU8R_2025-04-01Leuven_3-wave_longitudinal_study.csv) in folder temp, it will automatically use the OSF-hosted version to reproduce the results.

## Adherence to WORCS

This project uses the Workflow for Open Reproducible Code in Science (WORCS) to
ensure transparency and reproducibility. The workflow is designed to meet the
principles of Open Science throughout a research project. 

## More about WORCS

To learn how WORCS helps researchers meet the TOP-guidelines and FAIR principles,
read the preprint at https://osf.io/zcvbs/

# Miscellaneous

## Why is this repo named as the "islands of emotions"?

My training in Emotion-Focused Therapy (EFT) taught me that psychotherapy involves “islands of work in an ocean of empathy” [(Greenberg, 2017)](https://doi.org/10.1080/14779757.2017.1330702). This manuscript was inspired by that principle, and I am especially grateful to my EFT trainer, Ms. Elsie Lam (https://efthk.com/), as well as to the primary developers of EFT, including Les Greenberg. The central finding of this study—that getting through negative emotions involves transitioning between distinct negative emotions—can be seen as a reinterpretation of that guiding metaphor. To sail through the sea of negativity, one must first visit their islands of emotion.