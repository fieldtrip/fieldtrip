function y = negbin_pdf(x,l,r)
%NEGBIN_PDF Negative binomial probability density function
%
%  Description
%    Y = NEGBIN_PDF(X,L,R) returns the Negative binomial
%    probability density function with location parameter L and
%    dispersion parameter R (0<R<infty).
%
%    Negative binomial has different parameterizations and we use the form
%      p(x|l,r) = (r/(r+l))^r * gamma(r+y)/(gamma(r)*gamma(y+1))
%                             * (l/(r+l))^y
%    which approaches Poisson distribution when R goes to infinity.
%
%    The size of Y is the common size of X, L and R. A scalar input   
%    functions as a constant matrix of the same size as the other input.    
%
%    Note that the density function is zero unless X is an integer.
%
%    See also 
%      NEGBIN_LPDF, POISS_PDF, POISS_PDF

% Copyright (c) 2010 Jarno Vanhatalo, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

y=exp(negbin_lpdf(x,l,r));
