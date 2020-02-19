function test_issue1311

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_definetrial ft_preprocessing

%%

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.trialdef.eventtype = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue = [3 9]; % skip condition 5, to make it slightly more interesting
cfg = ft_definetrial(cfg);

trlfile = [tempname '.mat'];
trl     = cfg.trl;
save(trlfile, 'trl');


%%
% use the existing cfg

cfg.channel = 'STIM';
data0 = ft_preprocessing(cfg);


%%
% make a new one, with the trl from file

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.channel = 'STIM';

cfg.trl = trlfile;
data1 = ft_preprocessing(cfg);

% they should have the same data
assert(isequal(data0.time, data1.time));
assert(isequal(data0.trial, data1.trial));
% and the same sample and trial info
assert(isequal(data0.sampleinfo, data1.sampleinfo));
assert(isequal(data0.trialinfo, data1.trialinfo));
