# RegCombin: Partially Linear Regression under Data Combination.

We implement linear regression when the outcome of interest and some of the covariates are observed in two different datasets that cannot be linked, based on D'Haultfoeuille, Gaillac, Maurel (2022) <doi:10.3386/w29953>. 

The package allows for common regressors observed in both datasets, and for various shape constraints on the effect of covariates on the outcome of interest. It also provides the tools to perform a test of point identification. 

See the associated vignette https://github.com/cgaillac/RegCombin/blob/master/RegCombin_vignette.pdf for theory and code examples.

For a complete description of the method, see D’Haultfoeuille X, Gaillac C, Maurel A (2022). “Partially Linear Models under Data Combination” Working paper.
