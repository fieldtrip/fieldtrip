function v = logit(u)
%LOGIT Logit transformation
%
%  Description
%    V = LOGIT(U) computes the logit transformation of U
%
%      V = LOG(U./(1-U))
%  
%  See also LOGITINV
%
  
% Copyright (c) 2011 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

warning off MATLAB:divideByZero
v=reallog(u./(1-u));
warning on MATLAB:divideByZero
