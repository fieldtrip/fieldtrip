function [err, v] = gind2sub (sz, ndx)
%
%   [err,] = gind2sub ()
%
%Purpose:
%   This function converts a linear index to multiple subscripts based on the dimensions.
%   It is a further modified version from Ziad's modification based on the original Matlab 
%   function IND2SUB. The difference from what Ziad did is that we vary the subscript starting
%   from the last dimension and ends with the first dimension instead of the other way around.
%   
%Input Parameters:
%   sz: dimension vector
%   ndx: counter index for the total number of combinations.
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   v: array of subscripts  
%      
%Key Terms:
%   
%More Info :
%
%     Author : Gang Chen
%     Date : Tue Dec 23 10:57:48 EST 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda MD 20892


%Define the function name for easy referencing
FuncName = 'gind2sub';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

v = [ones(1,length(sz))];
if (length(ndx) ~= 1),
   fprintf(1,'Error zind2sub: length(IND) must be 1\n');
   v = [];
   return;
end

n = length(sz);
k = [1 cumprod(sz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1,
%for i = 1:1:n,
  v(i) = floor(ndx/k(i))+1;
  ndx = rem(ndx,k(i));
end

err = 0;
return;

