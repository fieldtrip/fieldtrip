function p=dir_pdf(x,a)
%DIR_PDF   Probability density function of uniform Dirichlet
%          distribution
%
%       Description:
%       P = DIR_PDF(X, A) returns the pdf of Dirichlet distribution 
%       with A at X

% Copyright (c) 2000 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


p=exp(gammaln(sum(a))-sum(gammaln)+sum(log(x).^(a-1)));
