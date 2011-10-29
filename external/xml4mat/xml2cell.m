function [y,varname]=xml2cell(filename_or_xml_string)

%XML2CELL reads non-MbML compliant xmlfile into matlab nested cell arrays
%
% Syntax [y,varname]=xmlfile2cell(filename_or_xml_string)
%
% Description:
%   1. Convert any non-MbML xml into MbML compliant string
%   2. Stores individual structures as nested cell arrays
%
% If it cannot be garanteed that non-MbML compliant XML has an
% internal referential consistency for convenient conversion to structures,
% this function builds the m-variable object model (MOM) as nested
% cell arrays. This approach ignores the possible dimensionality of the
% object and stores each entry, in a single cell, nested at the
% appropriate level.
% 
% Note 1 : if your XML string is MbML compliant use XML2MAT instead
% Note 2 : if your XML structure has consistent dimensionality use XML2STRUCT instead
%
% See also: XML2STRUCT, XML2MAT
%
% Jonas Almeida, almeidaj@musc.edu, 19 May 2003, MAT4NAT Tbox

% is this a xml string or xml filename ?
if filename_or_xml_string(1)=='<' 
    y=filename_or_xml_string;
else
    y=strrep(file2str(filename_or_xml_string),'''','''''');
end

% convert first to MbML compliant string and then onto an m-variable
[y,varname]=xml2mat(mbmling(y,0)); %ingnie: switched off displaying by making 0 of last argument