function [Sd] = zdeblank (S)
%
%   [Sd] = zdeblank (s)
%
%Purpose:
%   removes blanks surrounding a word in S
%
%
%Input Parameters:
%   S : a string
%
%
%Output Parameters:
%  Sd : S without the surrounding blanks
%
%
%Key Terms:
%
%More Info :
%   deblank, zdeblankall
%
%
%
%     Author : Ziad Saad
%     Date : Tue Sep 12 12:15:24 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'zdeblank';

blank_chars=sprintf(' \t\r\n\f\v%s',char(0));

non_blank_mask=~any(bsxfun(@eq,blank_chars,S(:)),2);
first_pos=find(non_blank_mask,1,'first');
last_pos=find(non_blank_mask,1,'last');

if isempty(first_pos)
    assert(isempty(last_pos));
    Sd='';
else
    Sd=S(first_pos:last_pos);
end
