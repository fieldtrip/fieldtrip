function test_bug3207

% WALLTIME 00:20:00
% MEM 12gb

% TEST ft_read_event read_edf

%% read data and annotations from a 2-channel test file

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3207/test3_2048Hz.EDF');

hdr1 = ft_read_header(filename, 'chanindx', 1);
dat1 = ft_read_data(filename, 'header', hdr1);
tim1 = (1:hdr1.nSamples)/hdr1.Fs;


hdr2 = ft_read_header(filename, 'chanindx', 2);
dat2 = ft_read_data(filename, 'header', hdr2);
tim2 = (1:hdr2.nSamples)/hdr2.Fs;

evt1 = ft_read_event(filename, 'detectflank', []); % expressed at 45 Hz
evt2 = evt1;
for i=1:numel(evt2)
  evt2(i).sample = evt2(i).timestamp*hdr2.Fs + 1;  % expressed at 2048 Hz
end

figure
hold on
plot(tim2, dat2 ./ max(dat2), 'g');
plot([evt1.timestamp], ones(1,numel(evt1)), 'ro');
plot([evt2.timestamp], ones(1,numel(evt2)), 'm+');

%% read the annotation events from the much larger file

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3207/ecog_mumc_anon.edf');

% use defaults
hdr = ft_read_header(filename);
dat = ft_read_data(filename);
tim = (1:hdr.nSamples)/hdr.Fs;

evt = ft_read_event(filename, 'detectflank', []);


%% do a quick inspection of the EDF data

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3207/ecog_mumc_anon.edf');
data = ft_preprocessing(cfg);

% read the events, don't detect flanks in a trigger channel but read annotations
event = ft_read_event(cfg.dataset, 'detectflank', []);

cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.event = event;
ft_databrowser(cfg, data);

%% do a full analysis on the EDF data

% this contains the data and the trialfun
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3207/'));

cfg = [];
cfg.dataset = 'ecog_mumc_anon.edf';
cfg.trialfun = 'trialfun_edf';
cfg = ft_definetrial(cfg);

cfg.demean = 'yes';
data = ft_preprocessing(cfg);

if false
  cfg = [];
  data = ft_rejectvisual(cfg, data);
end

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.layout = 'ordered';
ft_multiplotER(cfg, timelock);

cfg = [];
cfg.method = 'wavelet';
cfg.width = 5;
cfg.pad = 'nextpow2';
cfg.foilim = [1 70];
cfg.toi = -0.250:0.050:0.750;
cfg.trials = 1:10;
freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.layout = 'ordered';
cfg.baseline = [-inf 0];
cfg.baselinetype = 'relchange';
cfg.zlim = [0 5];
ft_multiplotTFR(cfg, freq);

