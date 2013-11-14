function y = geo_lpdf(x,p)
%GEO_LPDF     Geometric log probability density function (lpdf).
%   Y = GEO_LPDF(X,P) returns the log of geometric pdf with probability, P, 
%   at the values in X.
%
%   The size of Y is the common size of X and P. A scalar input   
%   functions as a constant matrix of the same size as the other input.    

% Copyright (c) 1999-2000 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

y = log(p) + log(1-p) * x;
