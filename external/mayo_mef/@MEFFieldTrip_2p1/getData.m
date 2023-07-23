function dat = getData(this, varargin)
% MEFFIELDTRIP_2P1.GETDATA read data from MEF 2.1 dataset for FieldTrip
%
% Syntax:
%   dat = getData(this)
%   dat = getData(__, channames, sampleunit, begsample, endsample, chanindx)
% 
% Input(s):
%   this            - [obj] MEFFieldTrip_2p1 object
%   channames       - [string] (opt) channel names (default = all channels)
%   sampleunit      - [char] (opt) unit of begsample and endsample (default
%                     = index)
%   begsample       - [num] (opt) time point of begin sample (default = 1st one)
%   endsample       - [num] (opt) time point of end sample (default = last one)
%   chanindx        - [num] (opt) index of selection from channames
%
% Output(s):
%   dat             - [num] data array channel x samples 
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2020 Richard J. Cui. Created: Thu 04/02/2020  4:38:22.701 PM
% $Revision: 0.1 $  $Date: Thu 04/02/2020  4:38:22.701 PM $
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
channames   = q.channames;
sampleunit  = q.sampleunit;
begsample   = q.begsample;
endsample   = q.endsample;
chanindx    = q.chanindx;

if isempty(channames)
    channames = this.ChannelName;
end % if

if isempty(begsample)
    begsample = 1;
end % if

if isempty(endsample)
    endsample = this.Samples;
end % if

if isempty(chanindx)
    chanindx = 1:numel(channames);
end % if

% =========================================================================
% main
% =========================================================================
sel_chan = channames(chanindx);
dat = this.importSession([begsample, endsample], sampleunit, this.SessionPath,...
    'SelectedChannel', sel_chan);

end % function getData

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% default
default_cn = "";
default_su = 'index';
default_bs = [];
default_es = [];
default_ci = [];

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('channames', default_cn, @isstring);
p.addOptional('sampleunit', default_su, @ischar);
p.addOptional('begsample', default_bs, @isnumeric);
p.addOptional('endsample', default_es, @isnumeric);
p.addOptional('chanindx', default_ci, @isnumeric);

% parse and return
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]
