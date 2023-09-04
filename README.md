# BvCategorical
An R package for fitting the bivariate categorical response regression model in high dimensions, as described in by [A likelihood-based approach for multivariate categorical response regression](https://www.tandfonline.com/doi/abs/10.1080/01621459.2021.1999819?journalCode=uasa20)]. Please contact [amolstad@umn.edu](mailto:amolstad@umn.edu) with any questions or comments. 

### Installation
BvCategorical can be loaded directly into R through the the `devtools` package:
```
install.packages("devtools")
library(devtools)
devtools::install_github("ajmolstad/BvCategorical")
library(BvCategorical)
```

### Citation instructions
If you use any of the functions from the package, please cite the article linked above. The bibtex entry is given below. 
```
@article{Molstad2023Likelihood,
  author = {Aaron J. Molstad and Adam J. Rothman},
  title = {A Likelihood-Based Approach for Multivariate Categorical Response Regression in High Dimensions},
  journal = {Journal of the American Statistical Association},
  volume = {118},
  number = {542},
  pages = {1402-1414},
  year  = {2023},
  publisher = {Taylor & Francis},
  doi = {10.1080/01621459.2021.1999819},
}
```
### Usage directions
Please visit [this example page]([https://ajmolstad.github.io/docs/HierMultinomExample.html](https://ajmolstad.github.io/docs/BvCategorical_Example.html)) for details on implementation and usage. 
