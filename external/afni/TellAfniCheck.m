function [err, g, b] = TellAfniCheck (w)
%
%   [err, g, b] = TellAfniCheck (w)
%
%Purpose:
%   Parses AFNI's returned string from plugout_drive
%   
%   
%Input Parameters:
%   w the string output by plugout_drive
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : N  N Problems
%   g   : Number of OK tokens
%   b   : Number of BAD tokens
%   
%      
%More Info :
%   
%    See TellAfni
%   
%
%     Author : Ziad Saad
%     Date : Wed Dec 7 11:00:06 EST 2005
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'TellAfniCheck';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
g = 0; 
b = 0;

slook = 'AFNI response string:';
nslook = length(slook);
[T] = strfind(w,slook);
nT = length(T);
if (isempty(T)),
   fprintf(2,'Warning: Response key not found in:\n%s\n', w);
   err = 0;
   return;
end

g = 0; b = 0; err = 0;
T(nT+1) = length(w);
for (i=1:1:nT),
   scons = w(T(i)+nslook:nslook+min(T(i)+4,T(i+1)));
   if (~isempty(strfind(scons,'OK'))) g = g + 1;
   elseif (~isempty(strfind(scons,'BAD'))) b = b + 1;
   else err = err + 1;
   end
end

return;

