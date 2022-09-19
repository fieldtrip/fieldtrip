function test_ft_write_data

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_write_data

%%

originaldir = pwd;
outputdir = tempname;
mkdir(outputdir);
cd(outputdir);

dataformat = {
  %  'fcdc_buffer'      % only works for online
  %  'fcdc_mysql'       % only works for online
  %   'mff'             % requires EEGLAB
  %   'neuralynx_ncs'   % requires FirstTimeStamp and TimeStampPerSample
  %   'homer_nirs'      % NIRS is tested further down
  %   'snirf'           % NIRS is tested further down
  'anywave_ades'
  'edf'
  'empty'
  'fcdc_global'
  'fcdc_matbin'
  'gdf'
  'matlab'
  'neuralynx_sdma'
  'plexon_nex'
  'brainvision_eeg'
  'brainvision_vhdr'
  'flac'
  'm4a' % not on linux
  'mp4' % not on linux
  'oga'
  'ogg'
  'wav'
  };

if ~ismac && ~ispc
  % not supported on linux
  dataformat = setdiff(dataformat, 'mp4');
  dataformat = setdiff(dataformat, 'm4a');
end


%%
% try the more general formats

fsample  = 44100; % this sampling rate is required for mp4 and m4a
nchan    = 1;
nsamples = 10*fsample;

dat = rand(nchan, nsamples); % between 0 and 1, otherwise the audio clips

hdr = [];
hdr.Fs          = fsample;
hdr.nChans      = nchan;
hdr.nSamples    = nsamples;
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;
for i=1:nchan
  hdr.label{i}       = num2str(i);
  hdr.chantype{i}    = 'unknown';
  hdr.chanunit{i}    = 'unknown';
end

for i=1:numel(dataformat)
  filename = fullfile(outputdir, ['write_' dataformat{i} '.ext']);
  ft_write_data(filename, dat, 'header', hdr, 'dataformat', dataformat{i});
end

%%
% try the NIRS formats

fsample  = 50;
nchan    = 10;
nsamples = 10*fsample;

dat = randn(nchan, nsamples);

opto = [];
opto.label         = {'Tx1-Rx1 680nm', 'Tx1-Rx1 690nm'};
opto.chanpos       = [0 0 1.5];
opto.optopos       = [0 0 1; 0 0 2; 0 0 3; 0 0 4];
opto.unit          = 'cm';
opto.optotype      = {'transmitter', 'receiver'};
opto.optolabel     = {'Tx1', 'Rx1', 'Tx2', 'Rx2'};
opto.wavelength    = [680 690];
opto.tra           = [
  +1 -1 0 0 % 1st and 2nd optode, 1st wavelength
  +2 -2 0 0 % 1st and 2nd optode, 2nd wavelength
  ];

hdr = [];
hdr.Fs          = fsample;
hdr.nChans      = nchan;
hdr.nSamples    = nsamples;
hdr.nSamplesPre = 0;
hdr.nTrials     = 1;
for i=1:nchan
  hdr.label{i}       = num2str(i);
  hdr.chantype{i}    = 'unknown';
  hdr.chanunit{i}    = 'unknown';
end

hdr.opto = opto;

hdr.chantype{1} = 'nirs';
hdr.chantype{2} = 'nirs';

% snirf
filename = fullfile(outputdir, 'write_snirf.snirf');
ft_write_data(filename, dat, 'header', hdr, 'dataformat', 'snirf');

% homer
filename = fullfile(outputdir, 'write_nirs.nirs');
ft_write_data(filename, dat, 'header', hdr, 'dataformat', 'homer_nirs');

%%
% clean up

cd(originaldir);
rmdir(outputdir, 's');
