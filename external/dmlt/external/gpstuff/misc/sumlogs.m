function y = sumlogs(x)
%SUMLOGS Sum of vector where numbers are represented by their logarithms.
%
%  Description
%    C=ADDLOGS(A) computes C=log(sum(exp(A))) in such a fashion
%    that it works even when elements have large magnitude.

% Copyright (c) 2003 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

y=x(1);
for k=2:length(x)
  y=addlogs(y,x(k));
end
