% TEST: ft_channelrepair fixsens test_ft_channelrepair

datainfo = test_datasets;

% get an MEG and an EEG set (hard-coded
eeginfo = datainfo(3);
meginfo = datainfo(7);

% do the MEG processing
fname = [meginfo.origdir,'raw/',meginfo.type,'preproc_',meginfo.datatype];
load(fname);

cfg = [];
cfg.method = 'triangulation';
neighbours = ft_neighbourselection(cfg, data);

cfg = [];
cfg.badchannel = data.label(100);
cfg.neighbours = neighbours;
newdata = ft_channelrepair(cfg, data);

% % do the EEG processing: this does not work, because there's no example EEG data
% with sensor positions 
% fname = [eeginfo.origdir,'raw/',eeginfo.type,'preproc_',eeginfo.datatype];
% load(fname);
% 
% cfg = [];
% cfg.method = 'template';
% cfg.template = 'EEG1010_neighb.mat';
% neighbours = ft_neighbourselection(cfg);
% 
% cfg = [];
% cfg.badchannel = data.label(10);
% cfg.neighbours = neighbours;
% newdata = ft_channelrepair(cfg, data);
% 
% 

