function c=addlogs(a,b)
%ADDLOGS Add numbers represented by their logarithms.
%
%  Description
%    C=ADDLOGS(A,B) computes C=log(exp(A)+exp(B)) in such a fashion
%    that it works even when A and B have large magnitude.

% Copyright (c) 2003 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if a>b
  c = a + log(1+exp(b-a));
else
  c = b + log(1+exp(a-b));
end
