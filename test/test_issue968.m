function test_issue968

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

p = tempname;
mkdir(p);
cd(p);

%%

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
cfg.trl = [1 30*300 0];
cfg.demean = 'yes';
cfg.baselinewindow = [-Inf 0];
cfg.channel = 'all';
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);

%%

hdr = ft_fetch_header(data);
dat = ft_fetch_data(data);

%%

filename = 'Subject01.ades';
ft_write_data(filename, dat, 'header', hdr, 'dataformat', 'anywave_ades');

%%
% the following is a fully simulated dataset
% it contains a block signal that is 0 at t=0, and goes up and down every second
% the amplitude is "realistic", i.e. 1e-6 Volt and 1e-12 Tesla

hdr = [];

hdr.label = {
  'EEG01'
  'ECG01'
  'EMG01'
  'SEEG01'
  'MEG01'
  'GRAD01'
  'Reference01'
  'Trigger01'
  };
hdr.chantype = {
  'eeg'
  'ecg'
  'emg'
  'seeg'
  'megmag'
  'megplanar'
  'refmag'
  'trigger'
  };
hdr.chanunit = {
  'V'
  'V'
  'V'
  'V'
  'T'
  'T/m'
  'T'
  'unknown'
  };
hdr.Fs = 300;

tim = (1:(30*300))/hdr.Fs;
dat = [
  mod(round(tim), 2) * 1e-6;
  mod(round(tim), 2) * 1e-6;
  mod(round(tim), 2) * 1e-6;
  mod(round(tim), 2) * 1e-6;
  mod(round(tim), 2) * 1e-12;
  mod(round(tim), 2) * 1e-12;
  mod(round(tim), 2) * 1e-12;
  mod(round(tim), 2);
  ];

filename = 'simulated.ades';
ft_write_data(filename, dat, 'header', hdr, 'dataformat', 'anywave_ades');

