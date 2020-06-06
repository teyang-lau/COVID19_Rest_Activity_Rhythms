# Rest Activity Rhythms During COVID19 Pandemic
Author: TeYang, Lau <br>
Last Updated: 6 June 2020

<img src = './Pictures/covid.jpg'>

<br>

### **Please refer to this [document]() for a more detailed description, analysis and insights of the project.** 
### **Alternatively, refer to the research journal article [Ong, Lau, Massar et al. 2020]() published in ... for the full scope and results of the project.

## **Project Goals** ##
1. To **investigate** rest activity rhythm changes throughout the COVID19 pandemic
2. **Cluster** days of steps counts to **identify** basis sets of rest activity rhythm (RAR) profiles
3. **Identify** groups of individuals who were differentially impacted by COVID19
4. **Identify** the sociodemographic composition of each group
5. **Investigate** how each group was affected in terms of sleep and physical activity 

## **Project Overview** ##
* Cleaning Fitbit data and filtering out abnormal days
* K-means++ clustering of 120,000+ days of intraday steps to identify sets of RAR profiles
* Computed proportion of time spent in each RAR profile per individual
* Hierarchical clustering of these proportions to identify groups of individuals differentially impacted by COVID19
* Performed Chi-Square Tests to look at sociodemographics in each group
* Performed Mixed ANOVAS to look at sleep and physical activity changes from before Circuit Breaker to during Circuit Breaker for each each group

## **About this dataset** ##
This data is shared by courtesy of the Singapore Health Promotion Board, as part of the [Health Insights Singapore (hiSG\)](https://www.hpb.gov.sg/hisg) study. hiSG is a population health study which seeks to better understand the health behaviours and lifestyles of Singapore residents through wearable technology. Participants were given devices (Fitbit™ Ionic, Fitbit™ Inc, San Francisco, CA) to track their activity/sleep and installed a mobile application to complete surveys over a period of 2 years. They were rewarded with points convertible to vouchers if they wore the tracker daily, logged sleep, meals, and completed surveys and were allowed to keep the device conditional on meeting study requirements.

hiSG consisted of a few phases targeting different population. The dataset used for this current analysis was from the first phase, which targeted young adults working in the Singapore Central Business District aged 21-40 years old.

Data was shared with the Sleep and Cognition Laboratory, Centre for Sleep and Cognition, Yong Loo Lin School of Medicine, National University of Singapore. This analysis is only a portion of the analysis that is conducted in the Sleep and Cognition Laboratory.





