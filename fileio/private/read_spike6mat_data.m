function dat = read_spike6mat_data(filename, varargin)

% read_spike6mat_data() - read Matlab files exported from Spike 6
%
% Usage:
%   >> header = read_spike6mat_data(filename, varargin);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'begsample'      first sample to read
%   'endsample'      last sample to read
%   'chanindx'  -    list with channel indices to read
%   'header'    - FILEIO structure header
%
% Outputs:
%   dat    - data over the specified range
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL
% Vladimir Litvak

if nargin < 1
    help read_spike6mat_data;
    return;
end

header    = ft_getopt(varargin, 'header');
begsample = ft_getopt(varargin, 'begsample');
endsample = ft_getopt(varargin, 'endsample');
chanindx  = ft_getopt(varargin, 'chanindx');

if isempty(header)
    header = read_spike6mat_header(filename);
end

if isempty(begsample), begsample = 1; end
if isempty(endsample), endsample = header.nSamples; end

try
    vars = struct2cell(load(filename));
catch
    ft_error('File not found or wrong format.');
end

if isempty(chanindx)
    chanindx = 1:numel(vars);
end

dat = zeros(length(chanindx), endsample-begsample+1);

for i = 1:length(chanindx)
    dat(i, :) = vars{chanindx(i)}.values(begsample:endsample);
end
