function test_fieldtrip2fiff()

if true % TODO: use file location on Donders server
  dataset_ctf = '/data1/projects/gosd/scripts/gosd/matlab/test_fieldtrip2fiff/Subject01/Subject01.ds';
  dataset_neuromag = '/data1/projects/gosd/scripts/gosd/matlab/test_fieldtrip2fiff/sample_audvis_raw.fif';
  dataset_eeg = '/data1/projects/gosd/recordings/svui/subjects/0003/eeg/raw/svui_0003_eeg_go-sd_010.raw';
  sensfile = '/data1/toolbox/fieldtrip/template/electrode/GSN-HydroCel-257.sfp';
  datadir = '/data1/projects/gosd/scripts/gosd/matlab/test_fieldtrip2fiff';
  
else 
  dataset_ctf = dccnfilename('/home/common/matlab/fieldtrip/data/Subject01.ds');
  dataset_neuromag = dccnfilename('?'); % MNE sample dataset
  dataset_eeg = dccnfilename('?'); % any EEG dataset
  datadir = dccnfilename('/home/common/matlab/fieldtrip/data/ftp/tutorial/test_fieldtrip2fiff'); % TODO: mkdir?
  
end

%-------------------------------------%
%-CTF: write raw
cfg = [];
cfg.dataset = dataset_ctf;
cfg.trialdef.triallength = Inf;
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';
cfg.channel = {'MEG', '-MLP31', '-MLO12'};
data = ft_preprocessing(cfg);

fifffile = [datadir filesep 'ctf_raw.fif'];
eventfile = [datadir filesep 'ctf_raw-eve.fif'];
fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
ft_read_header(fifffile);
ft_read_data(fifffile);
mne_read_events(eventfile); % this is simply a binary file with three columns
%-----------------%
%-------------------------------------%

%-------------------------------------%
%-CTF: write epoched
cfg = [];
cfg.dataset                 = dataset_ctf;   % name of CTF dataset  
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 3;                    % event value of FIC
cfg = ft_definetrial(cfg);            

cfg.channel   = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
data = ft_preprocessing(cfg);

fifffile = [datadir filesep 'ctf_epch.fif'];
fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
ft_read_header(fifffile);
ft_read_data(fifffile);
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
ft_read_header(fifffile);
ft_read_data(fifffile);
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
ft_read_header(fifffile);
ft_read_data(fifffile);
mne_read_events(eventfile); % this is simply a binary file with three columns
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
fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
ft_read_header(fifffile);
ft_read_data(fifffile);
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
ft_read_header(fifffile);
ft_read_data(fifffile);
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
hdr = ft_read_header(fifffile); 
ft_read_data(fifffile);
mne_read_events(eventfile); % this is simply a binary file with three columns
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
fieldtrip2fiff(fifffile, data)

%-----------------%
%-read back in
ft_read_header(fifffile);
ft_read_data(fifffile);
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
ft_read_header(fifffile);
ft_read_data(fifffile);
%-----------------%
%-------------------------------------%
