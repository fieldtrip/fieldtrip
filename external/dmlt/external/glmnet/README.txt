Glmnet for Matlab - Lasso (L1) and elastic-net regularized generalized linear models

Extremely efficient procedures for fitting the entire lasso or elastic-net regularization path for linear regression, logistic and multinomial regression models. The algorithm uses cyclical coordinate descent in a pathwise fashion 

Algorithm designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani. 
Fortran code written by Jerome Friedman
MATLAB wrapper written and maintained by Hui Jiang, jiangh@stanford.edu
Department of Statistics, Stanford University, Stanford, California, USA.

References

Regularized Paths for Generalized Linear Models via Coordinate Descent, Journal of Statistical Software, 33(1), 2010 
Pathwise coordinate optimization, The Annals of Applied Statistics. Volume 1, Number 2 (2007), 302-332.

=======================================================================

* Read MANUAL.txt for usage and examples.

* Right now only dense matrix is supported. We may add support for sparse matrix in the future.

* In this package, the following precompiled mex files are provided:

glmnetMex.mexw32, compiled with Windows XP and Matlab R2009b. When use, delete glmnetMex.dll to remove the warnings
glmnetMex.mexw64, compiled with Windows XP 64-bit, tested with Matlab R2009a and R2010a
glmnetMex.dll, compiled with Windows XP and Matlab R13, tested with Windows Vista and Matlab R14
glmnetMex.mexglx, compiled linux-ubuntu 9.10 x86-32bit
glmnetMex.mexa64, compiled with Linux x86-64bit and Matlab R2008b 64-bit
glmnetMex.mexmaci64, compiled with Mac OS X 10.5.8 Intel 64-bit and Matlab R2010a 64-bit (the beta trial version).

We may add precompiled versions for other platforms such as 32 bit Mac OS X in the future. For now, users of those platforms can easily compile the code themselves. To compile, simply type "mex glmnetMex.f glmnet.f" after you configured the fortran compiler in matlab, which depends on the compiler you are using. For example, you can use gfortran on mac OS X or g95 on Linux. For more details, you can refer to the matlab online help http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_external/f23674.html

* If you have successfully compiled a mex file which is not in the above list and would like to share it with other users, please contact us. We will be very appreciated for your help.

=======================================================================

Acknowledgement

glmnetMex.mexa64 was compiled with the help of Balasubramanian Narasimhan from Stanford University.
glmnetMex.mexmaci64 was kindly provided by Adam A. Margolin from the Broad Institute of Harvard and MIT.
glmnetMex.mexglx, cvglmnet.m and cvglmnetPlot.m were kindly provided by Bjorn Skovlund Dissing from Technical University of Denmark

Bugs reported by:
20 Oct 2009    Ramon Casanov from Wake Forest University
26 Jan 2010    Peter Rijnbeek from Erasmus University
09 Mar 2010    Ramon Casanov from Wake Forest University
25 Mar 2010    Gerald Quon from University of Toronto
