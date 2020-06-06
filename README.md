# Rest Activity Rhythms During COVID19 Pandemic
Author: TeYang, Lau <br>
Last Updated: 6 June 2020

<img src = './Pictures/covid.jpg'>

<br>

### **Please refer to this [document]() for a more detailed description, analysis and insights of the project.** ###

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
This data is shared by courtesy of the Singapore Health Promotion Board, as part of the [Health Insights Singapore (hiSG\)] (https://www.hpb.gov.sg/hisg) study. 
[Convert.ToBase64String(byte\[\])](http://msdn.microsoft.com/en-us/library/dhx0d524(v=vs.110).aspx)
Sleep and Cognition Laboratory, Centre for Sleep and Cognition, Yong Loo Lin School of Medicine, National University of Singapore

## **EDA** ##
Among other things, I looked at multicollinearity of the features and some of the relationship with age. Although correlation between variables are not very high, multicollinearity accessed by the variance inflation factor (VIF) is very high.

<img src = './Pictures/corr_hm.png' width='400'><img src = './Pictures/scatter_reg.png' width='400'>

## **Outlier removal using mahalanobis distance** ##
I chose to used mahalanobis distance to remove outliers as it is a multivariate distance measure and more suited for datasets with multiple features.

<img src = './Pictures/MD.png' width='400'>

## **Feature Importance** ##
I looked at the features that contribute most to predicting heart disease from the random forest and gradient boosting classifiers.
Thalassemia, number of major blood vessel, type of chest pain and ST depression are the top 4 features.

<img src = './Pictures/GB_featureImp.png' width='400'><img src = './Pictures/RD_featureImp.png' width='400'>

## **Model Comparisons** ##
Logistic Regression and Support Vector Machine performed the best on this dataset, both achieving an F1 score of 84.6%, precision of 88% and recall of 81.5%. Random Forest came in second with a F1 score of 83% while K-Nearest Neighbour performed the worst.

<img src = './Pictures/modelcompare.png'>

## **Conclusion** ##
To sum up the project, we attempted to predict heart disease using different machine learning models that differ in complexity. Out of all the models, logistic regression and support vector machines performed the best with F1 score of 85% and accuracy of 87%. In the era of big data, this sample size of 300 is considered quite small. With a much more richer dataset, I believe that a better model can be trained.

The second goal was to look at feature importance. Because of high multicollienarity in the features, the feature importance that we obtained from these models should be interpreted cautiously. Nevertheless, it appears that **thalassemia** is consistently the best predictor of heart disease. The number of major blood vessels that are blocked, the type of chest pains, and ST_depression are also good predictors of heart disease.

For a more detailed analyis and description of the project, please refer to this [notebook](https://nbviewer.jupyter.org/github/teyang-lau/Heart_Disease_Prediction/blob/master/predictors-of-heart-disease-an-exploration.ipynb)




