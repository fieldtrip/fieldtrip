function value=moxunit_fieldtrip_util_parse_dccnpath(line, unused)
% see if dccnpath is used as an expression
%
% value=moxunit_fieldtrip_util_parse_dccnpath(line)
%
% Input:
%   line                string which may contain a matlab expression
%                       involving dccnpath
%
% Output:
%   value               true if 'dccnpath' is present in line and not
%                       commented out or in a string, [] otherwise
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    if ~ischar(line)
        error('input must be numeric');
    end

    clean_line=line;

    % remove everything in quotes
    clean_line=regexprep(clean_line,'''.*''','');

    % remove everything after comment character
    clean_line=regexprep(clean_line,'%.*$','');

    % surround with spaces so we can match on non-word characters
    clean_line=[' ' clean_line ' '];

    % look for 'dccnpath' surrounded by non-word characters
    pat='\Wdccnpath\W';

    if isempty(regexp(clean_line,pat,'once'))
        value=[];
    else
        value=true;
    end



