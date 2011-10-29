function [y,varname]=xml2struct(filename_or_xml_string)

%XML2STRUCT reads non-MbML compliant xmlfile into matlab structure
%
% Syntax: [y,varname]=xml2struct(filename_or_xml_string)
%
% Description:
%   1. Convert any non-MbML xml into MbML compliant string
%   2. Stores consecutive structures in the a dimensional strucure
%
% If the non-MbML compliant XML has a consistent internal reference structure
% (those that were derived from explicit data models often do)
% this conversion will produce the best results, by building 
% 
% Note: if your XML string is MbML compliant use XML2MAT instead
%
% See also: XML2STRUCT
%
% Jonas Almeida, almeidaj@musc.edu, 19 May 2003, MAT4NAT Tbox

[y,varname]=xml2cell(filename_or_xml_string);
y=consolidateall(y);