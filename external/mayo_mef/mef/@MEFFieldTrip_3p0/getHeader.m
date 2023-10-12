function hdr = getHeader(this, channames)
% MEFFIELDTRIP_3P0.GETHEADER get header information of MEF 3.0 session
% 
% Syntax:
%   hdr = getHeader(this)
%   hdr = getHeader(this, channames)
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
% See also .

% Copyright 2020 Richard J. Cui. Created: Sun 03/22/2020  2:18:33.364 PM
% $Revision: 0.1 $  $Date: Sun 03/22/2020  2:18:33.364 PM $
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

% construct header structure
% --------------------------
orig = this.MetaData;

% sampling frequency
hdr.Fs = orig.time_series_metadata.section_2.sampling_frequency;

% number of channels
hdr.nChans = orig.number_of_time_series_channels;

% number of samples per trial
hdr.nSamples = orig.time_series_metadata.section_2.number_of_samples;

% number of trials
hdr.nTrials = 1;

% label of channel
l = convertStringsToChars(this.ChannelName);
hdr.label = l(:);

% channel type 
num_ch = hdr.nChans;
ct = cell(num_ch, 1);
ct(:) = {'unknown'};
hdr.chantype = ct;

% channel unit
unit_des = orig.time_series_metadata.section_2.units_description;
cu = cell(num_ch, 1);
cu(:) = {unit_des};
hdr.chanunit = cu;

% number of pre-trigger samples in each trial
hdr.nSamplesPre = this.Samples;

% save the original
hdr.orig = orig;

end % function

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