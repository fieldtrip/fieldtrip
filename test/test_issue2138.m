function test_issue2138

% WALLTIME 00:10:00
% MEM 2gb

%%

bidsroot = dccnpath('/home/common/matlab/fieldtrip/data/test/issue2138');

% this is the original
filename1 = fullfile(bidsroot, 'sub-001/ses-01/eeg/sub-001_ses-01_task-WorkingMemory_run-1_eeg.set');

% this is the one with enriched HED tags
filename2 = fullfile(bidsroot, 'sub-002/ses-01/eeg/sub-002_ses-01_task-WorkingMemory_run-1_eeg.set');

%%
% just to check that this works

hdr1 = ft_read_header(filename1, 'readbids', false);
event1 = ft_read_event(filename1, 'readbids', false);

hdr2 = ft_read_header(filename2, 'readbids', false);
event2 = ft_read_event(filename2, 'readbids', false);

%%
% default is readbids=true, which causes the sidecars to be read

hdr1 = ft_read_header(filename1);
event1 = ft_read_event(filename1);

hdr2 = ft_read_header(filename2);
event2 = ft_read_event(filename2);

%%
% show the events

cfg = [];
cfg.dataset = filename1;
cfg.trialdef.eventtype = '?';
ft_definetrial(cfg);

cfg = [];
cfg.dataset = filename2;
cfg.trialdef.eventtype = '?';
ft_definetrial(cfg);

%%
% this is some random ERP analysis using "normal" events

cfg = [];
cfg.dataset = filename1;
cfg.trialdef.eventtype = 'show_letter';
cfg.trialdef.eventvalue = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};
cfg.trialdef.prestim = 0.3;
cfg.trialdef.poststim = 0.7;
cfg = ft_definetrial(cfg);

cfg.baselinewindow = [-inf 0];
cfg.demean = 'yes';
rawdata = ft_preprocessing(cfg);

% make a selective average
cfg = [];
cfg.trials = ismember(rawdata.trialinfo.eventvalue, {'A' 'B' 'C'});
erp = ft_timelockanalysis(cfg, rawdata);

assert(sum(erp.cfg.trials)==9);

%%
% this is the same analysis for a dataset that also contains HED tags

cfg = [];
cfg.dataset = filename2;
cfg.trialdef.eventtype = 'show_letter';
cfg.trialdef.eventvalue = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};
cfg.trialdef.prestim = 0.3;
cfg.trialdef.poststim = 0.7;
cfg = ft_definetrial(cfg);

cfg.baselinewindow = [-inf 0];
cfg.demean = 'yes';
rawdata = ft_preprocessing(cfg);

% make a selective average
cfg = [];
cfg.trials = ismember(rawdata.trialinfo.eventvalue, {'A' 'B' 'C'});
erp = ft_timelockanalysis(cfg, rawdata);

assert(sum(erp.cfg.trials)==9);

%%
% this is a similar analysis, but now using HED tags

cfg = [];
cfg.dataset = filename2;
cfg.trialdef.eventtype = 'HED';
cfg.trialdef.eventvalue = 'Sensory-event,Visual-presentation,(Character,(Center-of,Computer-screen)),Target,*';
cfg.trialdef.prestim = 0.3;
cfg.trialdef.poststim = 0.7;
cfg = ft_definetrial(cfg);

cfg.baselinewindow = [-inf 0];
cfg.demean = 'yes';
rawdata = ft_preprocessing(cfg);

% make selective averages
cfg = [];
cfg.trials = contains(table2cell(rawdata.trialinfo), 'Items-to-memorize,Item-count,Target,Label/3');
erp_load3 = ft_timelockanalysis(cfg, rawdata);

cfg.trials = contains(table2cell(rawdata.trialinfo), 'Items-to-memorize,Item-count,Target,Label/5');
erp_load5 = ft_timelockanalysis(cfg, rawdata);

cfg.trials = contains(table2cell(rawdata.trialinfo), 'Items-to-memorize,Item-count,Target,Label/7');
erp_load7 = ft_timelockanalysis(cfg, rawdata);
