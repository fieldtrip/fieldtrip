function y = t_pdf(x, v, mu, sigma)
%T_PDF     Student's t probability density function (pdf).
%
%   Y = T_PDF(X,V,MU,SIGMA) Returns the Student's t pdf with
%   V degrees of freedom, MU mean and SIGMA standard deviation at X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

% Copyright (c) 1998-2008 Aki Vehtari

% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

if nargin < 4,
  sigma = 1;
end

if nargin < 3;
  mu = 0;
end

if nargin < 2,
  error('Requires at least two input arguments.');
end

term = gammaln((v + 1) / 2) - gammaln(v/2) -log(v.*pi)/2;
y = term - log(sigma) + log(1 + ((((x-mu)./sigma).^2)./v)) .* (-(v+1)/2);
y = exp(y);
