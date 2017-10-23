function header = read_spmeeg_header(filename)

% read_spmeeg_header() - import SPM5 and SPM8 meeg datasets
%
% Usage:
%   >> header = read_spmeeg_header(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   header   - FILEIO toolbox type structure
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Vladimir Litvak

if nargin < 1
  help read_spmeeg_header;
  return;
end;

try
    D = load(filename, 'D');
    D = D.D;
catch
    ft_error('File not found or wrong format.');
end

header = [];

if isfield(D, 'Radc') % SPM5
    header.Fs          = D.Radc;
    header.nChans      = D.Nchannels;
    header.nSamples    = D.Nsamples;
    
    if isfield(D, 'events') && isfield(D.events, 'start')
        header.nSamplesPre = D.events.start;
    else
        header.nSamplesPre = 0;
    end
    
    header.nTrials     = D.Nevents;
    header.label       = D.channels.name';
    
elseif all(isfield(D, {'type', 'Nsamples', 'Fsample', 'timeOnset'})) % SPM8

    header.Fs          = D.Fsample;
    header.nChans      = numel(D.channels);
    header.nSamples    = D.Nsamples;
    header.nSamplesPre = -D.timeOnset*D.Fsample;
    header.nTrials     = numel(D.trials);
    header.label       = {D.channels.label}';
    header.chantype    = lower({D.channels.type}');
    header.chanunit    = {D.channels.units}';
    
    if isfield(D.sensors, 'eeg') && ~isempty(D.sensors.eeg)
        header.elec = D.sensors.eeg;
    end
    
    if isfield(D.sensors, 'meg') && ~isempty(D.sensors.meg)
        header.grad = D.sensors.meg;
    end
    
else
    ft_error('Cannot recognize an SPM EEG header format');
end
    
header.orig = D;
