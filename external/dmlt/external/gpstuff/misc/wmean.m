function y = wmean(x, w)
% WMEAN   Weighted average or mean value.
%
%    Description
%      WMEAN(X,W) is the weighted mean value of the elements in X given 
%       weights W. 
%
%    See also wprctile, mean

% BUGS: Accepts only vector valued X

% Copyright (c) 2000-2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

y=sum(w.*x);
