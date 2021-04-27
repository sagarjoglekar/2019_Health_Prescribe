# 2019_Health_Prescribe
NHS prescription analysis, models and measurements

This directory is a living mess of things. There are notebooks which are written to test ideas, which were never used. 
Other notebooks were developed simply because someone thought it has to be done. 
All in all this project is going to scar you for a very long time. But just for the sake of maintaining some sanity in parsing the use of each notebook, Lets try and document the different "patterns" followed. 

## Requirements: 
1. The most important bit is the data. You would need to download the giant NHS prescriptions data files from [here](https://digital.nhs.uk/data-and-information/publications/statistical/practice-level-prescribing-data). There is no quick way around it. You will need to download each file my hand. There is no pattern in file urls, that could allow batch download. 

2. The dependencies are mostly vanilla, Pandas, networkx, statsmodel, geopandas etc. 

3. The paths to the downloaded data would need to be adjusted accordingly. 

I might have missed a whole lot of stuff, since for quick turn around I have used the magic store commands to pass objects between notebooks. I know, this is abominable sin and I have already made a deal with the devil to hire me as a software engineering manager in hell to torture programmer souls. 


## Approximate workflow. 

1. The most important script which does a bit of cleaning and parsing is 