Last modified: 2013-03-20 16:01:42 EET
-----------------------------------------------------------------

GPstuff: Gaussian process models for Bayesian analysis 4.0

Maintainers: Aki Vehtari <aki.vehtari@aalto.fi>
             Jarno Vanhatalo <jarno.vanhatalo@helsinki.fi>
             

This software is distributed under the GNU General Public Licence
(version 3 or later); please refer to the file Licence.txt,
included with the software, for details.

Table of contents:

1. INTRODUCTION
2. INSTALLING THE TOOLBOX
3. CONTENTS
4. TESTING THE INSTALLATION
5. USER GUIDE (VERY SHORT)
6. REFERENCE
7. JMLR DISCLAIMER

------------------------------------------
1. INTRODUCTION

  GPstuff is a collection of Matlab functions to build and analyze
  Bayesian models build over Gaussian processes. The toolbox is tested
  with Matlab R2009b-2012b in 64bit Windows and Linux environments and
  with Octave 3.6.4 in 64bit Linux environment (it should work in the
  32bit versions as well but they are not tested properly). NOTE! For
  Octave, use the Octave specific package available from the GPstuff
  homepage.

  The code for GPstuff can be found in subfolders. The SuiteSparse
  folder contains an exact copy of the SuiteSparse v3.4 toolbox by Tim
  Davis:
  http://www.cise.ufl.edu/research/sparse/SuiteSparse/current/SuiteSparse/
  The SuiteSparse is needed when using compactly supported covariance
  functions.

2. INSTALLING THE TOOLBOX

  If Matlab or Octave is started in root directory of GPstuff,
  startup.m script will add GPstuff subdirectories to the
  path. Alternatively, see startup.m for paths to add.

  Some of the functions in GPstuff are implemented using C in order to
  make the computations faster. In order to use these functions you
  need to compile them first. There are two ways to do that:

  1) Basic installation without compactly supported covariance
     functions

  * Install the GPstuff package by running matlab_install in this
    folder

  * With this option you are able to use all the other functions
    except for gpcf_ppcs*


  2) Installation with compactly supported covariance functions
  
  Compactly supported (CS) covariance functions are functions that
  produce sparse covariance matrices (matrices with zero elements). To
  use these functions (gpcf_ppcs*) you need the sparse GP
  functionalities in GPstuff which are build over SuiteSparse
  toolbox. To take full advantage of the CS covariance functions
  install GPstuff by running matlab_install('SuiteSparseOn' ) in the
  present directory.

    The function matlab_install compiles the mex-files and prints on
    the screen, which directories should be added to Matlab paths. 
    
3. CONTENTS
   
   The GPstuff packge contains the following subdirectories:
   diag  dist  gp  mc  misc  optim  xunit  SuiteSparse

   Each folder contains Contents.m, which summarizes the functions
   in the folder. 

   The 'gp' folder contains the main functionalities and demonstration
   programs. Other folders contain additional help functions.

4. TESTING THE INSTALLATION

   Installation can be tested by running command test_all, which
   runs all demos and compares the computed results to pre-saved
   results. Running test_all takes about one hour and it requires
   that 'xunit' toolbox is in the Matlab path. xunit package can be
   downloaded from
   http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework

5. USER QUIDE (VERY SHORT)

   It easiest to learn to use the package by running the demos. It is
   advisable to open the demo files in text editor and run them line
   by line. The demos are documented so that user can follow what
   happens on each line.

   The basic structure of the program is as follows. The program
   consist of separate blocks, which are:

      Gaussian process model structure (GP):
                      This is a structure that contains all the
                      model information (see GP_SET) and information
                      on, which inference scheme is used. 

                      GP structure contains covariance function
                      structures (GPCF_*) and likelihood structures
                      (LIK_*). 

      Covariance function structure (GPCF):
                      This is a structure that contains all of the
                      covariance function information (see
                      e.g. GPCF_SEXP). The structure contains the
                      hyperparameter values, pointers to nested
                      functions that are related to the covariance
                      function (e.g. function to evaluate covariance
                      matrix) and hyperprior structure.

      Likelihood structure:
                      This is a structure that contains all of the
                      likelihood function information (see
                      e.g. likelih_probit). The structure contains the
                      likelihood parameter values and pointers to
                      nested functions that are related to the
                      likelihood function (e.g. log likelihood and its
                      derivatives).

      Inference utilities:
                      Inference utilities consist of functions that
                      are needed to make the posterior inference and
                      predictions. These include, among others,
		        GP_OPTIM - Find MAP estimate for hyperparameters
                        GP_MC - Markov chain Monte Carlo sampling
                        GP_IA - Integration approximations
			GP_PRED - Predictions with Gaussian Process

6. REFERENCE

   If you use GPstuff, please use the reference:

   Jarno Vanhatalo, Jaakko Riihimäki, Jouni Hartikainen, Pasi Jylänki,
   Ville Tolvanen, Aki Vehtari (2013). GPstuff: A Toolbox for Bayesian
   Modeling with Gaussian Processes. In Journal of Machine Learning
   Research, accepted for publication.

7. JMLR DISCLAIMER

   THIS SOURCE CODE IS SUPPLIED \AS IS" WITHOUT WARRANTY OF ANY KIND, AND
   ITS AUTHOR AND THE JOURNAL OF MACHINE LEARNING RESEARCH (JMLR) AND
   JMLR'S PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL WARRANTIES,
   INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY
   AND FITNESS FOR A PARTICULAR PURPOSE, AND ANYWARRANTIES OR NON
   INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR
   USE OF THIS SOURCE CODE, AND NEITHER THE AUTHOR NOR JMLR, NOR JMLR'S
   PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR DAMAGES OF ANY KIND
   RESULTING FROM ITS USE. Without limiting the generality of the
   foregoing, neither the author, nor JMLR, nor JMLR's publishers and
   distributors, warrant that the Source Code will be error-free, will
   operate without interruption, or will meet the needs of the user.
