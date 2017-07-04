function header = read_spike6mat_header(filename)

% read_spike6mat_header() - read Matlab files exported from Spike 6
%
% Usage:
%   >> header = read_spike6mat_header(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   header   - FILEIO toolbox type structure
% _______________________________________________________________________
% Copyright (C) 2008 Institute of Neurology, UCL
% Vladimir Litvak

if nargin < 1
  help read_spike6mat_header;
  return;
end;

try
    vars = struct2cell(load(filename));
catch
    ft_error('File not found or wrong format.');
end

header = [];
header.nChans      = length(vars);
header.label = {};
header.orig = {};

fsample = [];
onsets = [];
lengths = [];
for i = 1:numel(vars)
    fsample(i) = round(1./vars{i}.interval);
    onsets(i)  = 1e-3*round(1e3*vars{i}.times(1));
    lengths(i) = vars{i}.length;
    header.label{i} = vars{i}.title;
    header.orig{i} = rmfield(vars{i}, {'values', 'times'});
end

if length(unique(fsample))>1 || length(unique(onsets))>1 || length(unique(lengths))>1
    ft_error('Only files with identical channel parameters are supported');
end

header.Fs          = unique(fsample);
    
header.nSamples    = unique(lengths);

header.nSamplesPre = -round(unique(onsets)*header.Fs);

header.nTrials     = 1;

header.label       = header.label(:);
    
