function test_issue1311

% WALLTIME 00:10:00
% MEM 6gb
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
% use the existing cfg as usual

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

%%
% process it as continuous

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.continuous = 'yes';
cfg.channel = 'STIM';
data2 = ft_preprocessing(cfg);

assert(numel(data2.trial)==1);


%%
% cut the continuous data in memory into trials

cfg = [];
cfg.trl = trlfile;
data3 = ft_redefinetrial(cfg, data2);

% they should have the same data
assert(isequal(data0.time, data3.time));
assert(isequal(data0.trial, data3.trial));
% and the same sample and trial info
assert(isequal(data0.sampleinfo, data3.sampleinfo));
assert(isequal(data0.trialinfo, data3.trialinfo));

%%
% shuffle and then check whether it returns the original order

ntrial  = numel(data3.trial);
shuffle = randperm(ntrial);
data3.trial       = data3.trial(shuffle);
data3.time        = data3.time(shuffle);
data3.trialinfo   = data3.trialinfo(shuffle,:);
data3.sampleinfo  = data3.sampleinfo(shuffle,:);

cfg = [];
cfg.trl = trlfile;
data4 = ft_redefinetrial(cfg, data3);

% they should have the same data
assert(isequal(data0.time, data4.time));
assert(isequal(data0.trial, data4.trial));
% and the same sample and trial info
assert(isequal(data0.sampleinfo, data4.sampleinfo));
assert(isequal(data0.trialinfo, data4.trialinfo));


