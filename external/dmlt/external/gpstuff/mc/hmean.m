function y = hmean(x)
% HMEAN   1/mean(1/input)
%
%   Description:
%   Y = HMEAN(X) evaluates y=1./mean(1./x)

% Copyright (c) 1998 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

y=1./mean(1./x);
