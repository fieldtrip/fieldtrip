function json = spm_jsonread(filename, opts)
% JSON (JavaScript Object Notation) parser - a compiled routine
% FORMAT json = spm_jsonread(filename, opts)
% filename - name of a JSON file or JSON string
% json     - JSON structure
% opts     - structure of optional parameters:
%              replacementStyle: string to control how non-alphanumeric
%                characters are replaced {'underscore','hex','delete','nop'}
%                [Default: 'underscore']
% 
% References:
%   JSON Standard: http://www.json.org/
%   JSMN C parser: http://zserge.com/jsmn.html
%   jsondecode: http://www.mathworks.com/help/matlab/ref/jsondecode.html
%__________________________________________________________________________
% Copyright (C) 2015-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_jsonread.m 7045 2017-03-17 10:41:12Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_jsonread.c not compiled - see Makefile')
