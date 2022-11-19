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
rawdata1 = ft_preprocessing(cfg);

% make a selective average
cfg = [];
cfg.trials = ismember(rawdata1.trialinfo.eventvalue, {'A' 'B' 'C'});
erp1 = ft_timelockanalysis(cfg, rawdata1);

assert(sum(erp1.cfg.trials)==9);

%%
% this is exactly the same analysis for a dataset that also contains HED tags
% but it does not use those HED tags

cfg = [];
cfg.dataset = filename2;
cfg.trialdef.eventtype = 'show_letter';
cfg.trialdef.eventvalue = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};
cfg.trialdef.prestim = 0.3;
cfg.trialdef.poststim = 0.7;
cfg = ft_definetrial(cfg);

cfg.baselinewindow = [-inf 0];
cfg.demean = 'yes';
rawdata2 = ft_preprocessing(cfg);

% make a selective average
cfg = [];
cfg.trials = ismember(rawdata2.trialinfo.eventvalue, {'A' 'B' 'C'});
erp2 = ft_timelockanalysis(cfg, rawdata2);

assert(sum(erp2.cfg.trials)==9);

%%
% this is the same analysis, now using the HED tags

cfg = [];
cfg.dataset = filename2;
cfg.trialdef.eventtype = 'HED';
cfg.trialdef.eventvalue = 'Sensory-event,Visual-presentation,(Character,(Center-of,Computer-screen)),Target,*';
cfg.trialdef.prestim = 0.3;
cfg.trialdef.poststim = 0.7;
cfg = ft_definetrial(cfg);

cfg.baselinewindow = [-inf 0];
cfg.demean = 'yes';
rawdata3 = ft_preprocessing(cfg);

% make a selective average
cfg = [];
cfg.trials = ~cellfun(@isempty, regexp(rawdata3.trialinfo.eventvalue, '(Character,Label/[ABC])'));
erp3 = ft_timelockanalysis(cfg, rawdata3);

assert(sum(erp3.cfg.trials)==9);

%%
% now using a trialfun that is specialized in HED tags
% it will assemble the HED tags when not present in the file

cfg = [];
cfg.dataset = filename2;
cfg.trialfun = 'ft_trialfun_hed';
cfg.trialdef.regexp = '.*,Target,\(Character,Label/[ABC]\),.*';
cfg.trialdef.prestim = 0.3;
cfg.trialdef.poststim = 0.7;
cfg = ft_definetrial(cfg);

cfg.baselinewindow = [-inf 0];
cfg.demean = 'yes';
rawdata4 = ft_preprocessing(cfg);

% make an average, it would be possible to make a selection using rawdata4.trialinfo
cfg = [];
erp4 = ft_timelockanalysis(cfg, rawdata4);

%% 
% just to demonstrate that reading+filtering and segmenting can also be done the other way around

cfg = [];
cfg.dataset = filename2;
cfg.hpfreq = 0.2;
cfg.hpfilter = 'yes';
rawdata5_continuous = ft_preprocessing(cfg);

% determine the selection of epochs
cfg = [];
cfg.dataset = filename2;
cfg.trialfun = 'ft_trialfun_hed';
cfg.trialdef.regexp = '.*,Target,\(Character,Label/[ABC]\),.*';
cfg.trialdef.prestim = 0.3;
cfg.trialdef.poststim = 0.7;
cfg = ft_definetrial(cfg);

% make the selection of epochs
rawdata5_epoched = ft_redefinetrial(cfg, rawdata5_continuous);

% make an average
cfg = [];
erp5 = ft_timelockanalysis(cfg, rawdata5_epoched);

%%
% all ERPs should be the same

assert(isequal(erp1.avg, erp2.avg));
assert(isequal(erp1.avg, erp3.avg));
assert(isequal(erp1.avg, erp4.avg));
% assert(isequal(erp1.avg, erp5.avg)); % this has been preprocessed differently

