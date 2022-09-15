function [err, snl, Nlines, strt, stp] = GetNextLine (s,n)
%
%   [err, snl, Nlines, strt, stp] = GetNextLine (s,n)
%
%Purpose:
%   reads the nth line of a string
%   lines are separated by a new line character, of course
%
%Input Parameters:
%   s
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%
%   snl the nth string
%   if you exceed the number of lines in the string, an empty string is returned
%   Nlines is the total number of lines in s
%   strt, stp are such that snl = s(strt:stp)
%Key Terms:
%
%More Info :
%   good for parins ls lists e.g:
%   [tm,v] = unix('ls -C1');
%   [err,sn] = GetNextLine(v,1)
%   [err,sn] = GetNextLine(v,4) etc
%
%
%     Author : Ziad Saad
%     Date : Mon Oct 11 13:15:13 CDT 1999


%Define the function name for easy referencing
FuncName = 'GetNextLine';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
strt = 0; stp = 0;

inl = find(char(s) == 10);

Nlines = length(inl);
if (n <= Nlines),
	if (n == 1),
			strt = 1;
		else
		strt = inl(n-1)+1;
	end
	stp = inl(n) -1;
	snl=s(strt:stp);
else
		snl = '';
end

err = 0;
return;

