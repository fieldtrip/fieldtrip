function [err, Vo, DS] = AFNI_OrientCode (V)
%
%   [err, Vo, DirSign] = AFNI_OrientCode (V)
%
%Purpose:
%   Changes the orientation code from number to letters and vice versa
%   
%   
%Input Parameters:
%   V, a 3x1 vector of numbers or characters
%      like [0 2 4] or 'RPI'
%   
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   Vo, a 3x1 vector of numbers or characters corresponding to the input
%   DirSign: is a 3x1 vector of +1 or -1 indicating whether the ith 
%       dimension of Vo is of the same direction (+1) or not (-1) 
%       with AFNI's colinear equivalent in the R A I convention.
%       
%   
%      
%Key Terms:
%   
%More Info :
%   AFNI .HEAD files
%   BrikInfo
%   [err, Vo, DirSign] = AFNI_OrientCode ([0 3 4]); Vo, DirSign
%   or
%   [err, Vo, DirSign] = AFNI_OrientCode ('SAR');Vo, DirSign
%
%     Author : Ziad Saad
%     Date : Tue Sep 5 19:09:52 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'AFNI_OrientCode';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
Vo = [];DS = [];

if (length(V) ~= 3),
	err = ErrEval(FuncName,'Err_Bad code length');
	return;
end

if (ischar(V)),
	Vo = [-1 -1 -1];
	for (i=1:1:3),
		switch V(i),
			case 'R'
				Vo(i) = 0; 
				DS(i) = 1;
			case 'L'
				Vo(i) = 1;
				DS(i) = -1;
			case 'P'
				Vo(i) = 2;
				DS(i) = -1;
			case 'A'
				Vo(i) = 3;
				DS(i) = 1;
			case 'I'
				Vo(i) = 4;
				DS(i) = 1;
			case 'S'
				Vo(i) = 5;
				DS(i) = -1;
			otherwise,
				err = ErrEval(FuncName,'Err_Cannot understand Orientation code');
				return;
		end
	end
else
	Vo = ['-' '-' '-'];
	for (i=1:1:3),
		switch V(i),
			case 0
				Vo(i) = 'R';
				DS(i) = 1;
			case 1
				Vo(i) = 'L';
				DS(i) = -1;
			case 2
				Vo(i) = 'P';
				DS(i) = -1;
			case 3
				Vo(i) = 'A';
				DS(i) = 1;
			case 4
				Vo(i) = 'I';
				DS(i) = 1;
			case 5
				Vo(i) = 'S';
				DS(i) = -1;
			otherwise,
				err = ErrEval(FuncName,'Err_Cannot understand Orientation code');
				return;
		end
	end
end




err = 0;
return;

