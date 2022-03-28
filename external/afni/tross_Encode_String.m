function [err,Sc] = tross_Encode_String (S)
%
%   [err,Se] = tross_Encode_String (S)
%
%Purpose:
%   encodes a string a la Tom Ross for use in AFNI's HEADER
%
%
%Input Parameters:
%   S
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   Se : the encoded version of S
%
%
%Key Terms:
%
%More Info :
%   README.attributes
%
%   To decode an encoded string (Se) you can run
%   eval(['fprintf(1,''' Se ''');'])
%
%     Author : Ziad Saad
%     Date : Wed Apr 11 09:28:10 PDT 2001
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'tross_Encode_String';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

%to find out the ascii number of special characters do
	%s = sprintf('\t'); double(s)
CR = 13;
LF = 10;
QUOTE = 34; %double quote, that makes C go nuts
TAB = 9;
BEL = 7;
VTAB = 11;
BS = 8;

Sd = double(S); %change S to ascii numbers
N_Sd = length(Sd);
Sc = '';
for (i=1:1:N_Sd),
	switch Sd(i),
		case CR,
			Sc = sprintf('%s\\r', Sc);
		case LF,
			Sc = sprintf('%s\\n', Sc);
		case QUOTE,
			Sc = sprintf('%s\\"', Sc);
		case TAB,
			Sc = sprintf('%s\\t', Sc);
		case BEL,
			Sc = sprintf('%s\\a', Sc);
		case VTAB,
			Sc = sprintf('%s\\v', Sc);
		case BS,
			Sc = sprintf('%s\\b', Sc);
		otherwise,
			Sc = sprintf('%s%s', Sc, char(Sd(i)));
	end
end


err = 0;
return;

