function [CnoHelp] = SkipMatlabHelp (C)
%
%   [CnoHelp] = SkipMatlabHelp (C)
%
%Purpose:
%   removes the first set of comments from a character array
%   that would normally consiture matlab's help text (see example)
%   
%   
%Input Parameters:
%   C : a character array
%   
%   
%Output Parameters:
%   CnoHelp : C without the help comments
%   
%   
%      
%More Info :
%   C = sprintf ('%%Hello Baby\n%%Couci Couci Cooo\nAbc = 100\n%%Waldo\njon=6')
%   [CnoHelp] = SkipMatlabHelp (C)
%   
%  PurgeComments
%	NextString
%	FindChar
%
%     Author : Ziad Saad
%     Date : Sat Mar 27 14:27:36 CST 1999 


%Define the function name for easy referencing
FuncName = 'SkipMatlabHelp';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
CnoHelp = [];


if (C(1) == '%'),
	cpos = 1;
	Loc = 1;
	while (~isempty(Loc) & Loc(1) == 1 & ~isempty(C)),
		%get the 1st line
			[err, Cnext, cend] = NextString (C, 'NewLine', 1);
		%get ridd of 1st line
			nC = length(C);
			C = C(cend+1:nC);
		%search for first %
			[err, Loc] = FindChar (C, '%');
	end
	CnoHelp = C;
else
	CnoHelp = C;
end

err = 0;
return;

