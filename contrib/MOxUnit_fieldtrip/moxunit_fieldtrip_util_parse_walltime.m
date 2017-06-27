function value=moxunit_fieldtrip_util_parse_walltime(line)
% parse WALLTIME property in fieldtrip test
%
% value=moxunit_fieldtrip_util_parse_walltime(line)
%
% Input:
%   line                string which may represent time; supported are:
%                       - 9             numeric input
%                       - '9'           char input
%                       - '01:12:13     char input; hours, minutes, seconds
%
% Output:
%   value               time in seconds, or the emty array [] if no value
%                       is reprsented in the input
%
% #   For MOxUnit_fieldtrip's copyright information and license terms,   #
% #   see the COPYING file distributed with MOxUnit_fieldtrip.           #

    if ~ischar(line)
        error('input must be numeric');
    end

    % try to parse time in HH:MM:SS format (hours:minutes:seconds)
    h_m_s_pat='^\s*(\d+):(\d+):(\d+)\s*';
    m=regexp(line,h_m_s_pat,'tokens','once');

    if ~isempty(m)
        secs=cellfun(@str2num,m);

        value=sum(secs(:) .* [3600; 60; 1]);
        return;
    end

    % try to parse a single non-negative integer
    pat='^\s*(\d+)\W+\s*';
    m=regexp([line ' '],pat,'tokens','once');

    if ~isempty(m)
        value=str2num(m{1});
        return;
    end


    value=[];