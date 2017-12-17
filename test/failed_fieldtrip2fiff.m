function failed_fieldtrip2fiff

% MEM 11gb
% WALLTIME 00:10:00

% TEST test_fieldtrip2fiff
% TEST fieldtrip2fiff ft_read_header ft_read_data ft_read_event

% use file location on Donders server
dataset_ctf      = dccnpath('/home/common/matlab/fieldtrip/data/MarkusBraille.ds');
sensfile         = dccnpath('/home/common/matlab/fieldtrip/template/electrode/GSN-HydroCel-257.sfp');
dataset_neuromag = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2036/sample_audvis_raw.fif');
dataset_eeg      = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2036/svui_0003_eeg_go-sd_010.raw');
datadir          = tempname;
mkdir(datadir);
  
%-------------------------------------%
%-CTF: write raw
cfg = [];
cfg.dataset = dataset_ctf;
cfg.trialdef.triallength = Inf;
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);

fifffile  = [datadir filesep 'ctf_raw.fif'];
eventfile = [datadir filesep 'ctf_raw-eve.fif'];
fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
cfg = [];
cfg.dataset = fifffile;
data1 = ft_preprocessing(cfg);

events = mne_read_events(eventfile); % this is simply a binary file with three columns
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-CTF: write epoched
cfg = [];
cfg.dataset             = dataset_ctf;   % name of CTF dataset  
cfg.trialdef.eventtype  = 'backpanel trigger';
cfg.trialdef.eventvalue = 4;
cfg.trialdef.prestim    = 0.4;
cfg.trialdef.poststim   = 0.6;

cfg = ft_definetrial(cfg);            

cfg.channel   = 'MEG';
%  data = ft_preprocessing(cfg); ERROR: not implemented

% fifffile = [datadir filesep 'ctf_epch.fif'];  % not clean yet, it writes evoked
% fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
% cfg = [];
% cfg.dataset = fifffile;
% data1 = ft_preprocessing(cfg);
% ft_datatype(data1)
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-CTF: write time-locked
cfg = [];
avg = ft_timelockanalysis(cfg, data);

fifffile = [datadir filesep 'ctf_avg.fif'];
fieldtrip2fiff(fifffile, avg)

%-----------------%
%-read back in
cfg = [];
cfg.dataset = fifffile;
data1 = ft_preprocessing(cfg);
avg1 = ft_timelockanalysis([], data1);
ft_datatype(avg1)  % returns 'timelock'
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-NEUROMAG: write raw
cfg = [];
cfg.dataset = dataset_neuromag;
cfg.trialdef.triallength = Inf;
cfg = ft_definetrial(cfg);

cfg.channel   = {'MEG', 'EEG', '-MEG 2443'  '-EEG 053'};
data = ft_preprocessing(cfg);

fifffile = [datadir filesep 'neuromag_raw.fif'];
eventfile = [datadir filesep 'neuromag_raw-eve.fif'];
fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
cfg = [];
cfg.dataset = fifffile;
data1 = ft_preprocessing(cfg);

events = mne_read_events(eventfile); % this is simply a binary file with three columns
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-NEUROMAG: write epoched
cfg = [];
cfg.dataset                 = dataset_neuromag;
cfg.trialdef.eventtype      = 'STI 014';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 1;
cfg = ft_definetrial(cfg); 

cfg.channel   = {'MEG', 'EEG', '-MEG 2443'  '-EEG 053'};
data = ft_preprocessing(cfg);

fifffile = [datadir filesep 'neuromag_epch.fif'];
%  fieldtrip2fiff(fifffile, data) ERROR: not implemented

%-----------------%
%-read back in
% cfg = [];
% cfg.dataset = fifffile;
% data1 = ft_preprocessing(cfg);
% ft_datatype(data1)
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-NEUROMAG: write time-locked
cfg = [];
avg = ft_timelockanalysis(cfg, data);

fifffile = [datadir filesep 'neuromag_avg.fif'];
fieldtrip2fiff(fifffile, avg)

%-----------------%
%-read back in
cfg = [];
cfg.dataset = fifffile;
data1 = ft_preprocessing(cfg);
avg1 = ft_timelockanalysis([], data1);
ft_datatype(avg1) 
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-EEG: write raw
cfg = [];
cfg.dataset = dataset_eeg;
cfg.trialdef.triallength = Inf;
cfg = ft_definetrial(cfg);

cfg.channel   = {'EEG'};
data = ft_preprocessing(cfg);
data.elec = ft_read_sens(sensfile);
data.elec.label{end} = 'E257'; %TODO: rename if necessary

fifffile = [datadir filesep 'eeg_raw.fif'];
eventfile = [datadir filesep 'eeg_raw-eve.fif'];
fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
% This will soft-crash if EEG has more than 60 electrodes because of mne2grad,
% l.193. fieldtrip2fiff stores the electrode locations in field
% "chs.eeg_loc" instead of "dig"
cfg = [];
cfg.dataset = fifffile;
data1 = ft_preprocessing(cfg);

events = mne_read_events(eventfile); % this is simply a binary file with three columns
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-EEG: write epoched
cfg = [];
cfg.dataset                 = dataset_eeg;
cfg.trialdef.eventtype      = 'trigger'; % TODO: adapt to EEG dataset
cfg.trialdef.prestim        = .5;
cfg.trialdef.poststim       = 1;
cfg.trialdef.eventvalue     = 'Targ'; % TODO: adapt to EEG dataset
cfg = ft_definetrial(cfg); 

cfg.channel   = {'EEG'};
data = ft_preprocessing(cfg);

fifffile = [datadir filesep 'eeg_epch.fif'];
%  fieldtrip2fiff(fifffile, data) ERROR: not implemented

%-----------------%
%-read back in
% cfg = [];
% cfg.dataset = fifffile;
% data1 = ft_preprocessing(cfg);
% ft_datatype(data1)
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-EEG: write time-locked
cfg = [];
avg = ft_timelockanalysis(cfg, data);

fifffile = [datadir filesep 'eeg_avg.fif'];
fieldtrip2fiff(fifffile, avg)

%-----------------%
%-read back in
cfg = [];
cfg.dataset = fifffile;
data1 = ft_preprocessing(cfg);
avg1 = ft_timelockanalysis([], data1);
ft_datatype(avg1) 
%-----------------%
%-------------------------------------%
