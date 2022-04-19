function [Sd] = zdeblankall (S)
%
%   [Sd] = zdeblankall (s)
%
%Purpose:
%   removes all blanks in a word in S
%   blanks are characters that return true
%   for the function isspace.
%   This includes space, tab, new line etc.
%
%Input Parameters:
%   S : a string
%
%
%Output Parameters:
%  Sd : S without any blanks
%
%
%Key Terms:
%
%More Info :
%   isspace, zdeblank
%
%   function is nothing but:
%   Sd = S(find(~isspace(S)));
%
%     Author : Ziad Saad
%     Date : Wed Nov 10 13:05:51 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


Sd = S(find(~isspace(S)));

return;

