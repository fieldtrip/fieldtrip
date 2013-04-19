function y = negbinztr_pdf(x,l,r)
%NEGBINZTR_PDF Zero trunc. negative binomial log probability density function
%
%  Description
%    Y = NEGBINZTR_PDF(X,L,R) returns the log of zero truncated
%    Negative binomial probability density function with location
%    parameter L and dispersion parameter R (0<R<infty).
%
%    Negative binomial has different parameterizations and we use the form
%      p(x|l,r) = (r/(r+l))^r * gamma(r+y)/(gamma(r)*gamma(y+1))
%                             * (l/(r+l))^y
%    which approaches Poisson distribution when R goes to infinity.
%    Zero truncated Negative Binomial has p(x==0|l,r)=0.
%
%    The size of Y is the common size of X, L and R. A scalar input
%    functions as a constant matrix of the same size as the other
%    input.
%
%     Note that the density function is zero unless X is an
%     integer.
%
%  See also 
%    NEGBINZTR_LPDF, NEGBIN_PDF, NEGBIN_LPDF

% Copyright (c) 2010 Jarno Vanhatalo, Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

y=exp(negbinztr_lpdf(x,l,r));
