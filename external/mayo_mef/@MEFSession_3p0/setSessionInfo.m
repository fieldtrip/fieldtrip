function metadata = setSessionInfo(this, varargin)
% MEFSESSION_3P0.SETSESSIONINFO set session information of MEFSession_3p0
%
% Syntax:
%   setSessionInfo(this, sesspath, password)
%   setSessionInfo(__, sort_channel)
% 
% Input(s):
%   this            - [obj] MEFSession_3p0 object
%   sesspath        - [char] session path of MEF 3.0 data
%   password        - [struct] MEF 3.0 password structure (see MEFSession_3p0
%                     for the details)
%   sort_channel    - [char] (opt) sort channel according to either 'alphabet' of
%                     the channel names or 'number' of the acquisiton
%                     channel number (default = 'alphabet')
%
%
% Output(s):
%   metadata        - [struct] MEF 3.0 MetaData (see
%                     read_mef_session_metadata_3p0)
%
% Example:
%
% Note:
%
% References:
%
% See also read_mef_session_metadata_3p0.

% Copyright 2020 Richard J. Cui. Created: Sat 03/21/2020 11:09:38.008 PM
% $Revision: 0.3 $  $Date: Fri 04/10/2020 12:05:51.930 AM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
sesspath = q.sesspath;
password = q.password;
sort_channel = q.sort_channel;

% =========================================================================
% main
% =========================================================================
this.SessionPath = sesspath; % set session path directory
this.Password = password; % set password
this.get_sess_parts;
this.get_sessinfo;

% the names in property ChannelName
metadata = this.read_mef_session_metadata_3p0;

if strcmpi(sort_channel, 'number')
    % sort the order of time_series_channels according to
    % acquisition_channel_number
    n_tsc = metadata.number_of_time_series_channels;
    acn = zeros(1, n_tsc);
    for k = 1:n_tsc
        acn(k) = metadata.time_series_channels(k).metadata.section_2.acquisition_channel_number;
    end % for
    [~, idx] = sort(acn); % from smallest to the largest
    metadata.time_series_channels = metadata.time_series_channels(idx);
end % if

% make sure that sequence of 'time_series_channels' is the same as that of
% ChannelName
meta_ch_name = convertCharsToStrings({metadata.time_series_channels.name});
this.ChannelName = meta_ch_name;
this.MetaData = metadata;

end % function setSessionInfo

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% default
default_sc = 'alphabet';

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('sesspath', @ischar);
p.addRequired('password', @isstruct);
p.addOptional('sort_channel', default_sc, @ischar)

% parse and return
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]
