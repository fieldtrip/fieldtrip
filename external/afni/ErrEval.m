function Decis = ErrEval (FuncName, ErrCode);
%
%	Decis = ErrEval (FuncName, ErrCode);
%
%This function displays an error message based on ErrCode and 
%returns 1 if it's a lethal error (Err_) , or 0 if it's 
%an error or a warning (Wrn_)that does not require aborting a function
%
%
%ErrCode is a string made out of two parts
%  The first part is a 4 letter string containing
%   Err_ or Wrn_  depending on wether it's an error or
%   a warning.
%  The second part is a string that indicates the type of
%   error encountered. 
%
%Defined values of ErrCode are :
%  No						%No problem
%  BadInpSize			%One of the input variables has a bad size	
%  BadInOut          %Bad combination of input and output parameters
%  FewInp				%Not enough inputs
%  ManyInp				%Too many inputs
%  FewOut				%Not enough outputs
%  ManyOut				%Too many outputs
%  SubFunc				%An error has occurred in a function within the script
%  BadOpt            %Bad Option
%  BadOptCombo       %Bad Option Combination
%  InpFieldMiss      %A field is missing from one of the input structures
%  OptFieldMiss      %A fiels is missing from the Options structures
%  EmptyInp          %Input data is empty 
%  FileExist         %File exists
%  FileNotExist      %File does not exist
% To pass a generic message, replace the ErrCode by any ther string
%
% Here's a typical example of the function usage
% if (ErrEval ('LoadSliceData','Err_Image Type Not Supported')), err = 1;	return;end
% 
%
%		Ziad Saad Thu Mar 26 11:43:51 CST 1998


if (eq_str(ErrCode,'No')), %in case you send a No only to 
	Decis = 0;
	return;
end

if (length(ErrCode) < 5),
	fprintf (2,'\nError in ErrEval, ErrCode %s is ambiguous.\nReturning a 1 decision.\n\n\a',ErrCode);
	Decis = 1;
	return;
end


tmp = ErrCode(1:3);

switch tmp
	case 'Err'
		Decis = 1;
	case 'Wrn'
		Decis = 0;
	otherwise
		fprintf (2,'\nError in ErrEval, prefix %s from %s is ambiguous.\nReturning a 1 decision.\n\n\a',tmp, ErrCode);
		Decis = 1;
		return;
end

tmp = ErrCode(5:length(ErrCode));

switch tmp
	case 'No'
		Decis = 0;		%even if it was Err_No
		return;
	case 'BadInpSize'
		s = 'Bad size of at least one of the input parameters.';
	case 'BadInOut',
		s = 'Bad combination of input and output parameters.';
	case 'FewInp'
		s = 'Too few inputs.';
	case 'ManyInp'
		s = 'Too many inputs.';
	case 'FewOut'
		s = 'Too few outputs';
	case 'ManyOut'
		s = 'Too many outputs';
	case 'SubFunc'
		s = 'Error (or Warning) in a function within this function (or script).';
	case 'BadOpt'
		s = 'Bad Option';
	case 'InpFieldMiss',
		s = 'A field is missing from one of the input structures.';
	case 'OptFieldMiss',
		s = 'A field is missing from the Options structure.';
	case 'BadOptCombo',
		s = 'Bad Combinations of Options.';
	case 'EmptyInp',
		s = 'One or all of input data is empty.';
	case 'FileExist',
		s = 'File specified exists.';
	case 'FileNotExist',
		s = 'File specified does not exist.';
	otherwise
		s = tmp;
end

%For the error message
if (Decis),
	fprintf (2,'\a\nError in %s : %s\n\n',FuncName,s);
else
	fprintf (2,'\nWarning from %s : %s\n\n',FuncName,s);
end

return;
