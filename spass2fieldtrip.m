function [lfp, spike, stm, bhv] = spass2fieldtrip(dirname, varargin)

% SPASS2FIELDTRIP reads data from a set of SPASS data files and converts
% the contents into data structures that FieldTrip understands. Note that
% dependent on the SPASS data it might be required to change some
% hard-coded parameters inside this function.
%
% Use as
%   [lfp, spike, stm, bhv] = spass2fieldtrip(dirname)
% Optionally you can specify the sample rate as key-value pairs
%  'fsample_ana' - default 1000
%  'fsample_swa' - default 32000
%
% The specified directory should contain the SPASS files, and the files should have
% the same name as the directory.
%
% The swa and sti input file are combined into the spike output structure.
% For the rest of the data it is trivial how the input and output relate.
%
% For example, if you specify
%   [lfp, spike, bhv, stm] = spass2fieldtrip('jeb012a02')
% then the following files should exist:
%   'jeb012a02/jeb012a02.ana'
%   'jeb012a02/jeb012a02.swa'
%   'jeb012a02/jeb012a02.spi'
%   'jeb012a02/jeb012a02.stm'
%   'jeb012a02/jeb012a02.bhv'
%
% Subsequently you can analyze the data in fieldtrip, or write the spike
% waveforms to a nex file for offline sorting using
%   ft_write_spike('jeb012a02_ch1.nex', spike, 'dataformat', 'plexon_nex', 'chanindx', 1)
%   ft_write_spike('jeb012a02_ch2.nex', spike, 'dataformat', 'plexon_nex', 'chanindx', 2)
%   ft_write_spike('jeb012a02_ch3.nex', spike, 'dataformat', 'plexon_nex', 'chanindx', 3)
%
% See also NUTMEG2FIELDTRIP, LORETA2FIELDTRIP, FIELDTRIP2SPSS

% Copyright (C) 2007, Robert Oostenveld
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance

fsample_ana = ft_getopt(varargin, 'fsample_ana', 1000);
fsample_swa = ft_getopt(varargin, 'fsample_swa', 32000);

anafile = fullfile(dirname, [dirname '.ana']);
swafile = fullfile(dirname, [dirname '.swa']);
spifile = fullfile(dirname, [dirname '.spi']);
stmfile = fullfile(dirname, [dirname '.stm']);
bhvfile = fullfile(dirname, [dirname '.bhv']);

if ~exist(anafile, 'file'), error('the file "%s" does not exist', anafile); end
if ~exist(swafile, 'file'), error('the file "%s" does not exist', swafile); end
if ~exist(spifile, 'file'), error('the file "%s" does not exist', spifile); end
if ~exist(stmfile, 'file'), error('the file "%s" does not exist', stmfile); end
if ~exist(bhvfile, 'file'), error('the file "%s" does not exist', bhvfile); end

% read the data
fprintf('reading %s\n', anafile); ana = read_labview_dtlg(anafile, 'int16');
fprintf('reading %s\n', swafile); swa = read_labview_dtlg(swafile, 'int16');
fprintf('reading %s\n', spifile); spi = read_labview_dtlg(spifile, 'int32');
fprintf('reading %s\n', stmfile); stm = read_labview_dtlg(stmfile, 'int32');
fprintf('reading %s\n', bhvfile); bhv = read_labview_dtlg(bhvfile, 'int32');

% determine the number of trials
ntrials = numel(bhv.data{1});

nchans = numel(ana.data)./ntrials;
ana.data = reshape(ana.data, [nchans, ntrials]);

% prepare the continuous LFP data
lfp         = [];
lfp.trial   = {};
lfp.time    = {};
lfp.label   = {};
lfp.fsample = fsample_ana;
for i=1:ntrials
  tmp = cell2mat(ana.data(:,i)')';
  lfp.trial{i} = tmp;
  lfp.time{i}  = ((1:size(tmp,2))-1)/fsample_ana;
  nsamples(i)  = length(lfp.time{i});
end
for i=1:size(ana.data,1)
  lfp.label{i,1} = sprintf('chan%d', i);
end

% the data is trial based, try to estimate the time between subsequent
% trials, or better: the time between subsequent stimuli
isi = 10^ceil(log10(max(nsamples)+1));

nchans = numel(swa.data)./ntrials;
swa.data = reshape(swa.data, [nchans, ntrials]);
spi.data = reshape(spi.data, [nchans, ntrials]);

% prepare the spike data
spike           = [];
spike.label     = {};
spike.waveform  = {}; % 1xnchans cell-array, each element contains a matrix (Nsamples X Nspikes), can be empty
spike.timestamp = {}; % 1xnchans cell-array, each element contains a vector (1 X Nspikes)
spike.unit      = {}; % 1xnchans cell-array, each element contains a vector (1 X Nspikes)
for i=1:size(swa.data,1)
  spike.label{i} = sprintf('chan%d', i);
end

nspikes = zeros(nchans, ntrials);
for i=1:nchans
  for j=1:ntrials
    nspikes(i,j) = length(swa.data{i,j});
  end
end

% reinsert the inter-trial intervals, this links both the LFP and the spike
% timestamps to a common continuous timeaxis
for j=1:nchans
  for i=1:ntrials
    spi.data{j,i} = spi.data{j,i} + (i-1)*(fsample_swa./fsample_ana)*isi;
  end
end

% convert the spike timestamps and waveforms to a fieldtrip-compatible format
for i=1:nchans
  spike.waveform{i}  = cell2mat(swa.data(i,:));
  spike.timestamp{i} = cell2mat(spi.data(i,:)')';
  spike.unit{i}      = ones(size(spike.timestamp{i}));  % since unsorted
end

% prepare the other output
stm = stm.data{1}(:);
bhv = bhv.data{1}(:);

% store some additional information in the cfg structure
cfg = [];

% remember where the lfp trials are on the imaginary continuous timeaxis,
% this links both the LFP and the spike timestamps to a common continuous
% timeaxis
trl = zeros(ntrials, 3);
for i=1:ntrials
  begsample = (i-1)*isi+1;
  endsample = begsample+nsamples(i)-1;
  offset    = 0;
  trl(i,:) = [begsample endsample offset];
end
cfg.trl = trl;

% store the header information
lfp.hdr.FirstTimeStamp = 0;
lfp.hdr.TimeStampPerSample = fsample_swa./fsample_ana;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance lfp spike
ft_postamble history lfp spike
