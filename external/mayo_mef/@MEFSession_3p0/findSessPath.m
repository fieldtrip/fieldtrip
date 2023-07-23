function [sesspath, channames] = findSessPath(this, filename)
% MEFFSESSION_3P0.FINDSESSPATH find session path and channel name 
%
% Syntax:
%   [sesspath, chan_sel] = findSessPath(filename)
%
% Input(s):
%   this            - [obj] MEFFieldTrip_3p0 object
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

% Copyright 2020 Richard J. Cui. Created: Sun 03/22/2020  7:51:37.906 AM
% $Revision: 0.2 $  $Date: Thu 04/02/2020 10:41:14.141 AM $
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
[~, ~, f_ext] = fileparts(filename);
switch f_ext
    case '.mefd' % session folder
        if isfolder(filename)
            % session path
            sesspath = filename;
            % channel names
            listing = dir(fullfile(filename, '*.timd'));
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
        else
            error('MEFSession_3p0:FindSessPath:invalidFiletype',...
                'invalide file type. can only process MEF 3.0 files')
        end % if
    case '.timd' % channel folder
        if isfolder(filename)
            [sesspath, channames] = fileparts(filename);
        else
            error('MEFSession_3p0:FindSessPath:invalidFiletype',...
                'invalide file type. can only process MEF 3.0 files')
        end % if
    case '.segd' % segment folder
        if isfolder(filename)
            chanpath = fileparts(filename); % channel path
            [sesspath, channames] = fileparts(chanpath);
        else
            error('MEFSession_3p0:FindSessPath:invalidFiletype',...
                'invalide file type. can only process MEF 3.0 files')
        end % if
    case {'.tdat', '.tidx', '.tmet'} % data file
        if isfile(filename)
            segpath = fileparts(filename); % segment path
            chanpath = fileparts(segpath); % channel path
            [sesspath, channames] = fileparts(chanpath);
        else
            error('MEFSession_3p0:FindSessPath:invalidFiletype',...
                'invalide file type. can only process MEF 3.0 files')
        end % if
    otherwise
        error('MEFSession_3p0:FindSessPath:invalidFiletype',...
            'invalide file type. can only process MEF 3.0 files')
end % end switch
channames = convertCharsToStrings(channames);

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
