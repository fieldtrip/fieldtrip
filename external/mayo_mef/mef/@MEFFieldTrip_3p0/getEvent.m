function evt = getEvent(this, channames)
% MEFFIELDTRIP_3P0.GETEVENT get MEF 3.0 events for FieldTrip
%
% Syntax:
%   evt = getEvent(this)
%   evt = getEvent(this, channames)
% 
% Input(s)
%   this            - [obj] MEFFieldTrip_3p0 object
%   channames       - [string] channel names
% 
% Output(s)
%   hdr             - [struct] structure of header information (see
%                     ft_read_header for the detail)
% 
% Note:
% 
% References:
%
% See also .

% Copyright 2020 Richard J. Cui. Created: Sun 03/22/2020  3:10:07.514 PM
% $Revision: 0.1 $  $Date: Sun 03/22/2020  3:10:07.519 PM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, channames);
channames = q.channames;

% =========================================================================
% main
% =========================================================================
% check channel names if provided
% -------------------------------
if ~isempty(channames)
    sess_chan = this.ChannelName; % all channels in the session
    if all(ismember(channames, sess_chan)) == false
        error('MEFFieldTrip_3p0:getHeader:invalidChannel',...
            'invalid channel names')
    end % if
end % if

% construct the events
% ---------------------
cont = table2struct(this.SessionContinuity); % table of continuity

num_evt = numel(cont); % number of events
evt = struct([]);
for k = 1:num_evt
    % type
    evt(k).type = 'continuity recording';
    % sample
    evt(k).sample = cont(k).SampleIndexStart;
    % value
    evt(k).value = [];
    % offset
    evt(k).offset = cont(k).SampleIndexEnd;
    % duration
    evt(k).duration = cont(k).SegmentLength;
    % time stamp of sample start
    evt(k).timestamp = cont(k).SampleTimeStart;
end % for

end % function getEvent

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% default
default_cn = "";

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('channames', default_cn, @isstring);

% parse and return
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]
