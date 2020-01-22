# dsmatch
The R package *dsmatch* provides matching algorithm based on both propensity score and 
prognostic score to estimate average treatment effect and quantile 
treatment effect. Classical matching algortihms such as propensity score
matching and prognositc score matching are also contained for comparison.

## Installation with `devtools`
```R
devtools::install_github("Yunshu7/dsmatch")
```

## Main paper 
**Double score matching estimators of average and quantile treatment effects** available on arxiv 
<https://arxiv.org/abs/2001.06049>

## Main functions
There are four functions in the package:
- dsmatchATE: Double Score Matching Estimator for Average Treatment Effect
- dsmatchQTE: Double Score Matching Estimator for Quantile Treatment Effect
- dsmatchATT: Double Score Matching Estimator for Average Treatment Effect for the Treated
- dsmatchQTT: Double Score Matching Estimator for Quantile Treatment Effect for the Treated  

Detailed information and examples are included in the documents of each function. For example:
```R
?dsmatchATE
```










