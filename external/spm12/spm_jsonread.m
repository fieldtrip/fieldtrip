function json = spm_jsonread(filename, opts)
% JSON (JavaScript Object Notation) parser - a compiled routine
% FORMAT json = spm_jsonread(filename, opts)
% filename - name of a JSON file or JSON string
% json     - JSON structure
% opts     - structure or list of name/value pairs of optional parameters:
%              ReplacementStyle: string to control how non-alphanumeric
%                characters are replaced {'underscore','hex','delete','nop'}
%                [Default: 'underscore']
%              Prefix: string to prepend when first character of a field is
%                not alphabetical [Default: 'x']
% 
% References:
%   JSON Standard: https://www.json.org/
%   JSMN C parser: https://zserge.com/jsmn/
%   jsondecode: https://www.mathworks.com/help/matlab/ref/jsondecode.html
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2015-2023 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
%error('spm_jsonread.c not compiled - see Makefile')

persistent runonce
if nargin > 1 && isempty(runonce)
    warning('Optional parameters disabled when using builtin jsondecode.');
    runonce = 1;
end

json = jsondecode(fileread(filename));
