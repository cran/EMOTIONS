# EMOTIONS: Ensemble Models fOr lacTatION curveS

**INDEX**

-   📝 [Introduction](#introduction)
-   🛠️ [Installation](#Installing-the-Package)
-   💻 [Analysis](#Analysis)
-   🤝 [Authors](#Authors)
-   📜 [Licence](#License)
-   ✉️ [Contact](#Contact)

## **Introduction**

The mathematical representation of lactation curves has significant applications in various areas of animal science. Models that simulate milk production under different conditions are valuable for physiologists, nutritionists, and geneticists, allowing them to study mammary gland function and test hypotheses. These models also support management decisions related to timing and efficiency. Over the past decades, numerous models have been proposed for representing lactation curves in dairy species. These models differ mainly in their regression type (linear or nonlinear), number of parameters, relationships among parameters, and ability to represent lactation patterns such as peak yield, time at peak, and persistency.

Despite the advantages of having a diverse range of models, selecting a single model to represent the lactation curve can present limitations. Typically, model selection is based on comparing metrics such as the Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC). However, this metric-based approach assumes that the chosen metric is an optimal criterion for model selection. In practice, these metrics may fail to correctly identify the best model, especially when no single model is clearly superior. Furthermore, using a single model selected based on these metrics can lead to overfitting and introduce bias, particularly when the dataset is noisy or contains numerous variables.

Ensemble modeling and model averaging are powerful techniques that enhance robustness, accuracy, and generalization in predictive modeling. Instead of relying on a single model, ensemble methods combine multiple models to reduce variance, mitigate overfitting, and improve predictive performance. Model averaging incorporates model uncertainty by weighting predictions according to their posterior probabilities, leading to more reliable estimates. These approaches are particularly valuable when individual models exhibit varying performance across different datasets or conditions. By integrating multiple models, ensemble methods improve stability and resilience, making them especially useful in complex biological and ecological systems where uncertainty quantification is crucial.

The **EMOTIONS** package provides a set of tools for fitting 47 different lactation curve models previously reported in the literature. Some of these models and the pre-defined starting parameters were obtained from the lactcurves R package (https://CRAN.R-project.org/package=lactcurves). Once the data is fitted to each model, ensemble predictions are generated using bagging based on AIC, BIC, root mean square percentage error (RMSPE), mean absolute error (MAE), and variance. Additionally, the package provides predictions for daily milk records using Bayesian Model Averaging (BMA) and calculates cosine similarity for each model's predictions. The ranking of models across individual predictions can be visualized using the **RidgeModels** and **ModelRankRange** functions, which help users better understand the weight assigned to each model. Furthermore, the **PlotWeightLac** function allows users to compare predicted and actual daily milk records. Lastly, **EMOTIONS** enables the estimation of resilience indicators based on lag-1 autocorrelation, logarithm of residual variance, and residual skewness using predicted daily milking records.

## **Installing the Package**

Use the following code to install **EMOTIONS**:

```{r, eval=FALSE}
devtools::install_github("https://github.com/pablobio/EMOTIONS")
```

## **Analysis**

### **Loading the Package**

```{r}
library(EMOTIONS)
```
The vignette with the tutorial for the use of the main functions from EMOTIONS is available here: https://rpubs.com/pablo_bio/EMOTIONS_vignette

## Authors

- **Pablo A. S. Fonseca**  
  Instituto de Ganadería de Montaña (CSIC-Univ. de León)  
  24346 Grulleros, León  
  Email: [p.fonseca@csic.es](mailto:p.fonseca@csic.es)

- **Marcos Prates**  
  Department of Statistics, Universidade Federal de Minas Gerais, Brasil  

- **Aroa Suárez-Vega**  
  Dpto. Producción Animal, Facultad de Veterinaria, Universidad de León (24007)  

- **Ruth Arribas Gonzalo**  
  Dpto. Producción Animal, Facultad de Veterinaria, Universidad de León (24007)  

- **Beatriz Gutierrez-Gil**  
  Dpto. Producción Animal, Facultad de Veterinaria, Universidad de León (24007)  

- **Juan José Arranz**  
  Dpto. Producción Animal, Facultad de Veterinaria, Universidad de León (24007) 

## License

GPL-3


## Contact

**Pablo Fonseca**

[![Gmail Badge](https://img.shields.io/badge/-psouf@unileon.es-c14438?style=flat-square&logo=Gmail&logoColor=white&link=mailto:psouf@unileon.es)](mailto:psouf@unileon.es)
[![Google Scholar Badge](https://img.shields.io/badge/Google-Scholar-lightgrey)](https://scholar.google.com/citations?user=1VUm8EIAAAAJ&hl=pt-BR)
[![ResearchGate Badge](https://img.shields.io/badge/Research-Gate-9cf)](https://www.researchgate.net/profile/Pablo_Fonseca2)

<!-- display the social media buttons in your README -->


[![alt text][6.1]][6]


<!-- links to social media icons -->
<!-- no need to change these -->

<!-- icons with padding -->

[6.1]: http://i.imgur.com/0o48UoR.png (github icon with padding)

<!-- icons without padding -->

[6.2]: http://i.imgur.com/9I6NRUm.png (github icon without padding)


<!-- links to your social media accounts -->
<!-- update these accordingly -->

[6]: https://github.com/pablobio


<!-- Please don't remove this: Grab your social icons from https://github.com/carlsednaoui/gitsocial -->

**Grupo MEjoraGeneticaAnimal-ULE (MEGA_ULE)**

[![alt text][1.1]][1]

[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)

[1.2]: http://i.imgur.com/wWzX9uB.png (twitter icon without padding)

[1]: https://twitter.com/MEGA_ULE
