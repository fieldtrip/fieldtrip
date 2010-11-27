function [err, Cnext, cend] = NextString (C, Delim, cpos, Opt)
%
%   [err,Cnext, cend] = NextString (C, Delim, cpos, Opt)
% or [Cnext, cend] = NextString (C, Delim, cpos, Opt)
%Purpose:
%   gets the next string in charachter array C that is delimited by Delim
%   
%   
%Input Parameters:
%   C : a character array
%   Delim : delimitation character, see FindChar for special characters
%      Note that Delim is included in Cnext
%   cpos : position to start from in C ( cpos >=1 )
%  Opt optional options structure with the following fields
%    .Deblank (0/1) default 0, removes trailing white spaces (tab, newline, formfeed).
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%	Cnext : C(pos:'location of next Delim')   
%  cend is such that Cnext = C(cpos:cend);
%   Warning : cend is not affected for Deblanking, the equality on
%      the previous line is before deblanking is done on Cnext
%      
%More Info :
%   FindChar
%   SkipMatlabHelp
%   PurgeComments
%
% example
% ks = sprintf ('I am hyper, Yes Oh Yes I am.\nHow about You\tLaddy. Are You hyper ?\n')
% [err, Cnext, cend] = NextString (ks,'.',1)
% or 
% [err, Cnext, cend] = NextString (ks,'.',29)
% [err, Cnext, cend] = NextString (ks,'?',29)
% or 
% [err, Cnext, cend] = NextString (ks,'v',1)
% or
% [err, Cnext, cend] = NextString (ks,',',1)
% [err, Cnext, cend] = NextString (ks,',',12)
% or
% [err, Cnext, cend] = NextString (ks,'NewLine',1)
% [err, Cnext, cend] = NextString (ks,'Space',1)
%
% You get the idea right ?
%
%     Author : Ziad Saad
%     Date : Sat Mar 27 14:38:45 CST 1999 


%Define the function name for easy referencing
FuncName = 'NextString';

%Debug Flag
DBG = 1;

if (nargin == 3),
	Opt.Deblank = 0;
end

if (~isfield(Opt,'Deblank') | isempty(Opt.Deblank)), Opt.Deblank = 0; end

%initailize return variables
err = 1;
Cnext = '';

nC = length(C);

if (cpos > nC),
	err = ErrEval(FuncName,'Err_cpos > length(C)');
	return;
end

[err, Loc] = FindChar (C(cpos:nC), Delim);

if (isempty(Loc)),
	cend = nC-cpos;
else
	switch Delim
		case 'NewLine',
			offset = -1;
		case '\n',
			offset = -1;
		case 'Tab',
			offset = -1;	
		case '\t',
			offset = -1;
		case 'Space',
			offset = -1;
		case ' ',
			offset = -1;
		otherwise,
			offset = -1;
	end
	cend = Loc(1)+cpos+ offset;
end

Cnext = C(cpos:cend);
nCnext = length(Cnext);

if (Opt.Deblank & nCnext),
	ispc = isspace(Cnext);
	if (sum(ispc) == length(ispc)),
		%all space
		Cnext = '';
	else
		%get ridd of first blanks
		cnt = 1;
		while (ispc(cnt)),
			cnt = cnt + 1;
		end
		cnt2 = nCnext;
		while (ispc(cnt2)),
			cnt2 = cnt2 - 1;
		end
		Cnext = Cnext(cnt:cnt2);
	end
end

err = 0;

if (nargout < 3), 
   err = Cnext; 
   Cnext = cend; 
end

return;



