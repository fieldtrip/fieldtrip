function value=moxunit_fieldtrip_util_parse_dependency(line)
% parse WALLTIME property in fieldtrip test
%
% value=moxunit_fieldtrip_util_parse_walltime(line)
%
% Input:
%   line                string which represents one or more alphanumeric
%                       strings
%
% Output:
%   value               cell string with the strings in line separated by
%                       white space
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    if ~ischar(line)
        error('input must be numeric');
    end

    value=regexp(line,'(\w+)','match');

