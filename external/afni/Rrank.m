function p = Rrank(R,tol)
%
%   [err,] = Rrank.m ()
%
%Purpose:
%   
%   
%   
%Input Parameters:
%   
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   
%   
%      
%Key Terms:
%   
%More Info :
%   
%   
%   
%
%     Author : Gang Chen
%     Date : Tue Mar 23 16:31:43 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda MD 20892


%Define the function name for easy referencing
FuncName = 'Rrank.m';

%Debug Flag
DBG = 1;


if nargin<2
   tol = 100 * eps * max(size(R));
end
if (min(size(R))==1)
   d = abs(R(1,1));
else
   d = abs(diag(R));
end
p = sum(d > tol*max(d));

return;

