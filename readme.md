External Validation of a Mobile Clinical Decision Support System for Diarrhea Etiology Prediction
in Children: A Multicenter Study in Bangladesh and Mali

This prospective observational study aimed to develop and externally validate the accuracy of a mobile software application (“App”) for the prediction of viral-only etiology of acute diarrhea in children 0-59 months in Bangladesh and Mali using a previously derived and internally validated model from Brintz, Ben J., et al. "A modular approach to integrating multiple data sources into real-time clinical prediction for pediatric diarrhea." Elife 10 (2021): e63009.

Please contact Ben Brintz at ben.brintz@hsc.utah.edu for information regarding analyses. All analyses were conducted in R version 3.6.2. 

Annotated code (Phase_1_elife.R) and files necessary to complete the analysis are provided at https://github.com/LeungLab/DiaPR_Phase1. 

All human subjects are de-identified and use only a study id. 


Files of interest: 
DiaPRData2020.01.01.csv and base_smartphone.csv contain the CRF information from both Bangladesh and Mali. 
dlbg with AFEs.csv and Mali_Afe.csv contain the TAC data needed to define etiology
dhaka_weather.csv and Mali_weather.csv contain weather data needed for weather component of model 
Files with Mod(s), OR, and xs contain fitted model and odds ratio information needed to calculate post-test odds 