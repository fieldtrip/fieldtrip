function yml=simplify_mbml(xml)

%SIMPLIFY_MBML simplifies the MbML representation by removing attributes with default values
%Syntax: yml=simplify_mbml(xml)
%Description:
%   - Input arguments - 
%     xml is a MbML XML statement (string)
%   - Output arguments - 
%     yml is the simplified version (string)
%
%Jonas Almeida, almeidaj@musc.edu,29 Oct 2002

%1. remove size for char and double when they are null, elements or
%horizontal vectors

yml=regexprep(xml,' class=(("double")|("char")|("struct")|("cell")) size="[01] \d*"',' class=$1','tokenize');
yml=regexprep(yml,' class="char">','>');

