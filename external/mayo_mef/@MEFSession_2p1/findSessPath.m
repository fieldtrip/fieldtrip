function [sesspath, channames] = findSessPath(this, filename)
% MEFFSESSION_2P1.FINDSESSPATH find session path and channel name 
%
% Syntax:
%   [sesspath, chan_sel] = findSessPath(filename)
%
% Input(s):
%   this            - [obj] MEFFieldTrip_2p1 object
%   filename        - [char] name of the file or folder of the dataset
%
% Output(s):
%   sesspath        - [char] session path
%   channames       - [string] channel names included
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2020 Richard J. Cui. Created: Thu 04/02/2020 12:23:41.042 PM
% $Revision: 0.1 $  $Date: Thu 04/02/2020 12:23:41.042 PM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, filename);
filename = q.filename;

% =========================================================================
% main
% =========================================================================
if isempty(filename)
    sesspath = '';
    channames = "";
    return
end % if

% find session path and channel name
% ----------------------------------
if isfolder(filename)
    % get channel names
    listing = dir(fullfile(filename, '*.mef'));
    if isempty(listing)
        channames = ""; % empty string
    else
        num_chan = numel(listing);
        channames = strings(1, num_chan);
        for k = 1:num_chan
            chan_k = convertCharsToStrings(listing(k).name);
            [~, channame_k] = fileparts(chan_k);
            channames(k) = channame_k;
        end % for
    end % if
    sesspath = filename;
elseif isfile(filename)
    [sesspath, fname, f_ext] = fileparts(filename);
    if isequal(f_ext, '.mef')
        channames = convertCharsToStrings(fname);
    else
        channames ="";
    end % if
else
    error('MEFSession_2p1:FindSessPath:invalidFiletype',...
        'invalide file filder or type. can only process MEF 2.1 files')
end % if

end % function findSessPath

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% default

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('filename', @ischar);

% parse and return
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]
