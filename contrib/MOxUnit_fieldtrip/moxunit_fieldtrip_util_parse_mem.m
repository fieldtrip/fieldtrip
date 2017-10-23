function value=moxunit_fieldtrip_util_parse_mem(line)
% parse MEM property in FieldTrip test
%
% value=moxunit_fieldtrip_util_parse_mem(line)
%
% Input:
%   line                string which may representin integer. Suffixes of
%                       'kb', 'mb', 'gb' and 'tb' are suported
%                       - '9'           char input
%                       - '9kb'         kilobytes
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

    % try to parse
    pat='\s+(\d+|inf)\s*([tTgGmMkK][bB])?\W+';
    m=regexp([' ' line ' '],pat,'tokens','once');

    if ~isempty(m)
        value=str2num(m{1});

        if numel(m)<=1 || isempty(m{2})
            multiplier=1;
        else
            % kilo, mega, giga or tera
            powers='kmgt';
            i=find(lower(m{2}(1))==powers);
            assert(numel(i)==1);
            multiplier=2^(i*10);
        end

        value=value*multiplier;
        return;
    end

    value=[];
