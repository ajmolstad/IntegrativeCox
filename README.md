# IntegrativeCox
A repository for code accompanying the article [Dimension reduction for integrative survival analysis](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13736).  

### Contents
In the **Functions** directory, you will find **LRCox.R**, which contains the function for fitting the model proposed in the manuscript. Specifically, function **IntCox** fits the entire solution path for a vector of candidate *s* and vector of candidate *r*. The function **IntCoxCV** performs *K* fold cross-validation (and does not fit the solution path to the entire dataset). To use these functions, you must also source the C++ file **updateBeta.cpp**. The third script in the Functions directory, **RRCox_PPG.R** fits the nuclear norm and group lasso penalized estimator discussed in Appendix A using a proximal-proximal gradient descent algorithm. 

In the Simulations directory, one can find the scripts needed to reproduce exactly the results from Section 5 of the manuscript. Specifically, ModelX_Main.R creates data for one of three models, then sources the Fit_Main.R script, which fits the models and computes the performance metrics. Each bash script is used to initialize all replicates for each of the three models.  Of course, one will need to carefully modify file paths since this is taken directly from A. J. Molstad's working directory on HiperGator and the University of Florida. 

For clarifications or specific usage instructions, please contact [amolstad@umn.edu](amolstad@umn.edu). 

### Citation instructions
Please cite the most recent version of the article mentioned above. As of June 2023, this was the following (in bibtex): 
```
@article{Molstad2023Dimension,
  author = {Molstad, Aaron J. and Patra, Rohit K.},
  title = {Dimension reduction for integrative survival analysis},
  journal = {Biometrics},
  volume = {n/a},
  number = {n/a},
  pages = {},
  keywords = {Cox proportional hazards model, dimension reduction, integrative survival analysis, majorize-minimize, penalty method, reduced-rank regression, variable selection},
  doi = {https://doi.org/10.1111/biom.13736},
  url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13736},
  eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1111/biom.13736},
}
```
