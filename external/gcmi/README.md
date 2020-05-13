# gcmi : Gaussian-Copula Mutual Information

Functions for calculating mutual information and other information theoretic quantities using a parametric Gaussian copula. 

If you use this code for analysis that is published in an indexed journal or repository, please cite the following article:

RAA Ince, BL Giordano, C Kayser, GA Rousselet, J Gross and PG Schyns  
"A statistical framework for neuroimaging data analysis based on mutual information estimated via a Gaussian copula"  
*Human Brain Mapping* (2017) 38 p. 1541-1573 [doi:10.1002/hbm.23471](http://dx.doi.org/10.1002/hbm.23471)

For journals with supplementary information that may not be indexed for citations, **please place the citation in the indexed main manuscript**.

The [`matlab_examples`](matlab_examples) directory contains tutorial example scripts reproducing the analyses from that paper.

## Installation

#### Matlab

Add the contents of the `matlab` directory to your Matlab path.

#### Python

Put `gcmi.py` somewhere on your `PYTHONPATH` or in your working directory. `import gcmi` and then prefix the functions below with the module name, e.g. `gcmi.gcmi_cc(x,y)`.

## Usage

This is an overview of the functions available. Fuller descriptions of function arguments are provided in the in-line documentaion (Python docstrings, Matlab help strings).

The function names follow a format where the first part of the name is the quantity which the function calculates (eg mi, cmi, gcmi) followed by an underscore, followed by the types of variables that function acts on, in order corresponding to the arguments - g: Gaussian, d: Discrete, c: Continuous (any distribution). 

Please note a crucial difference in the interface between the Python and Matlab versions for multidimensional variables. In keeping with the convention of each platform in Matlab samples should be the first axis of the array (data stored in Fortran order); in Python samples should be the last axis of the array (data stored in C order by default).

| Platform | Multivariate data layout    |
| :------: | :----------------------:    |
| Matlab   | samples FIRST axis (rows)   |
| Python   | samples LAST axis (columns) |

Discrete inputs are passed as two arguments, the vector of values over samples y, and an integer parameter Ym specifying the size of the discrete space. In Python discrete variables should be stored in an integer data type array; in Matlab a double is used but should contain only integer values. y takes values between `0` and `Ym-1` inclusive. Empty classes are not supported.

Care should be taken with continuous variables that contain many repeated values. The copula transform which depends on a rank ordering will not be well defined in such cases. Possible approaches include repeated calculations while jittering the data with low amplitude noise to avoid the numerically equivalent values, or using binning and discrete methods. 

For functions with a `biascorrect` option, this is an optional true or false switch (default true) which indicates whether analytic bias correction for the entropy of Gaussian variables is applied. The bias correction increases computation time and is not needed when combined with permutation testing.

#### GCMI functions

`gcmi` functions estimate the Gaussian Copula Mutual Information, including input data checking and the copula transform step. We suggest new who are trying out the measure start with these functions.

*  `I = gcmi_cc(x,y)` 

    Calculate GCMI between two (possibly multidimensional) continuous variables x and y. x and y can have any marginal distribution but should not contain repeated values. 

*  `I = gcmi_model_cd(x,y,Ym)` 

    Calculate GCMI between a (possibly multidimensional) continuous variable x and a discrete y (with values between 0 and Ym-1 inclusive). This uses ANOVA style model comparison. It provides a lower bound to MI in the 1d case (but not in higher dimensions).

*  `I = gcmi_mixture_cd(x,y,Ym)` 

    Calculate GCMI between a (possibly multidimensional) continuous variable x and a discrete y (with values between 0 and Ym-1 inclusive). This uses a Gaussian mixture model, fit with robust measures of location and scale.

*  `I = gccmi_ccc(x,y,z)` 

    Calculate GCMI between two (possibly multidimensional) continuous variables x and y conditioned on a third, z. x / y / z can have any marginal distribution but should not contain repeated values.

*  `I = gccmi_ccd(x,y,z,Zm)` 

    Calculate GCMI between two (possibly multidimensional) continuous variables x and y conditioned on a third discrete variable, z (with values between 0 and Zm-1 inclusive). x and y can have any marginal distribution but should not contain repeated values.

#### Low level functions

These functions implement the different steps for the GCMI calculation. They are provided separately for computational efficiency (e.g. copula transform only needs to be performed once prior to permutation testing).

##### Copula transformation functions

*  `c = ctransform(x)`

    Compute empirical CDF value (copula transform) for each value. If x is >2D transformation is performed on each dimension separately. 

*  `cx = copnorm(x)` 

    Perform copula normalisation (equivalent to `norminv(ctransform(x))`). Returns standard normal samples with rank ordering preserved. If x is >2D normalization is performed on each dimension separately.  

##### Information theoretic quantities for Gaussian variables

These functions calculate information theoretic quantities (mi: mutual information, cmi: conditional mutual information) for Gaussian variables. Together with copula normalization above they implement the GCMI estimator: `gcmi(x,y) = mi_gg(copnorm(x),copnorm(y))`.

*  `I = mi_gg(x,y,biascorrect)` 

    Calculate MI between two (possibly multidimensional) Gaussian variables x and y. 

*  `I = mi_model_gd(x,y,Ym,biascorrect)` 

    Calculate MI between a (possibly multidimensional) Gaussian variable x and a discrete y (with values between 0 and Ym-1 inclusive). This uses ANOVA style model comparison. It provides a lower bound to MI in the 1d case (but not in higher dimensions). 

  *  `I = mi_mixture_gd(x,y,Ym,biascorrect)` 

    Calculate MI between a (possibly multidimensional) Gaussian variable x and a discrete y (with values between 0 and Ym-1 inclusive). This uses a Gaussian mixture model.

*  `I = cmi_ggg(x,y,biascorrect)` 

    Calculate CMI between two (possibly multidimensional) Gaussian variables x and y conditioned on a third, z. 

##### Miscellaneous functions

*  `H = ent_g(x, biascorrect)`

    Return analytic entropy of a (possibly multidimensional) Gaussian variable. 


## Compiled parallel versions

The code here implements the method in a direct, portable way with error checking.
If you find the measures work for you and plan an analysis of a large data set you might find the performance of this implementation a bottleneck.
I have a fast parallellised compiled mex implementation of the functions which can speed up typical neuroimaging mass-univariate analysis by two orders of magnitude, as well as access to excellent computing resources at the University of Glasgow.
If this sounds like it would be useful I would be very happy to hear from anyone interested in a collaboration.

## FAQ

### I get low or zero estimated MI values when there is a non-linear relationship in my data

GCMI is an approximate lower bound estimate to the true non-linear mutual information. GCMI quantifies dependence that can be modelled with a Gaussian copula, therefore can be thought of as having similar sensitivity to a rank correlation (in the univariate case). If you have a distribution with non-linear effects (e.g. a donut shape) the GCMI lower bound might not be tight and it could even give a value close to 0. If you need sensitivity to all forms of non-linear you should use a different estimator. The advantages of GCMI are its robustness, computational efficiency and performance for multivariate calculations, but it is a lower bound approximation. 

### I tried GCMI on my data but I get negative values. I thought MI was always positive?

We cannot measure the "true" value of the mutual information for the system studied, we can only estimate the quantity from the data we observe. Any estimate from data has a variance - if you repeat the experiment you will not get exactly the same value. For a positive quantity like MI this means that even when there is no dependence, you will get a different positive value for each data set. This means the estimator has a bias - the expected value when there is no effect (e.g. mean of permutations) is greater than zero. When we apply bias correction techniques to remove this bias and get the mean closer to zero, this means there will be more negative values observed when the variables are independent. Also, since the `mi_model_gd` method compares different models, which each feature different sampling variability, this estimator can have negative values even before bias correction.

![Figure: `mi_model_gd` null permutation histograms](matlab_examples/bias_demo.png)

The figure above shows the results of `mi_model_gd` calculated with 100 trials each of two classes, when the data in each class are sampled from the same standard normal distribution. This is repeated 1000 times (see [`bias_demo.m`](matlab_examples/bias_demo.m)). Vertical red lines show the mean of the repetitions. You can see that many simulations result in negative values. 

## Version history

* 0.3 - Change `mi_mixture` from an unconditional copnorm to a more faithful mixture model based on robust measures  of location and  scale
* 0.2 - Rename `mi_gd` functions to `mi_model_` and add `mi_mixture`
* 0.1 - Original version
