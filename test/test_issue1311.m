function test_issue1311

% WALLTIME 00:20:00
% MEM 6gb
% DEPENDENCY ft_definetrial ft_preprocessing

%%

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
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
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
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
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
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

%%
% process the artifacts with the data from disk

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
cfg.continuous = 'yes';
cfg.trl = trlfile;

cfg1 = cfg;
cfg1.artfctdef.clip.channel = 'MZC01';
[cfg1, artifact1] = ft_artifact_clip(cfg1);

cfg2 = cfg;
cfg2.artfctdef.ecg.channel = 'MZC01';
cfg2.artfctdef.ecg.feedback = 'no'; % otherwise it is interactive
[cfg2, artifact2] = ft_artifact_ecg(cfg2);

cfg3 = cfg;
cfg3.artfctdef.eog.channel = {'MLF*', 'MRF*'};
[cfg3, artifact3] = ft_artifact_eog(cfg3);

cfg4 = cfg;
cfg4.artfctdef.jump.channel = 'MZC01';
[cfg4, artifact4] = ft_artifact_jump(cfg4);

cfg5 = cfg;
cfg5.artfctdef.muscle.channel = {'MLT*', 'MRT*'};
[cfg5, artifact5] = ft_artifact_muscle(cfg5);

cfg6 = cfg;
cfg6.artfctdef.nan.channel = 'MZC01';
[cfg6, artifact6] = ft_artifact_nan(cfg6);

cfg7 = cfg;
cfg7.artfctdef.threshold.channel = 'MZC01';
[cfg7, artifact7] = ft_artifact_threshold(cfg7);

cfg8 = cfg;
cfg8.method = 'detect';
cfg8.artfctdef.tms.channel = 'MZC01';
[cfg8, artifact8] = ft_artifact_tms(cfg8);

cfg9 = cfg;
cfg9.artfctdef.zvalue.cutoff = 3;
cfg9.artfctdef.zvalue.channel = 'MZC01';
[cfg9, artifact9] = ft_artifact_zvalue(cfg9);

%%
% process the artifacts with the data from memory

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
cfg.continuous = 'yes';
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);

%%
cfg = [];
cfg.trl = trlfile;

cfg1 = cfg;
cfg1.artfctdef.clip.channel = 'MZC01';
[cfg1, artifact1] = ft_artifact_clip(cfg1, data);

cfg2 = cfg;
cfg2.artfctdef.ecg.channel = 'MZC01';
cfg2.artfctdef.ecg.feedback = 'no'; % otherwise it is interactive
[cfg2, artifact2] = ft_artifact_ecg(cfg2, data);

cfg3 = cfg;
cfg3.artfctdef.eog.channel = {'MLF*', 'MRF*'};
[cfg3, artifact3] = ft_artifact_eog(cfg3, data);

cfg4 = cfg;
cfg4.artfctdef.jump.channel = 'MZC01';
[cfg4, artifact4] = ft_artifact_jump(cfg4, data);

cfg5 = cfg;
cfg5.artfctdef.muscle.channel = {'MLT*', 'MRT*'};
[cfg5, artifact5] = ft_artifact_muscle(cfg5, data);

cfg6 = cfg;
cfg6.artfctdef.nan.channel = 'MZC01';
[cfg6, artifact6] = ft_artifact_nan(cfg6, data);

cfg7 = cfg;
cfg7.artfctdef.threshold.channel = 'MZC01';
[cfg7, artifact7] = ft_artifact_threshold(cfg7, data);

cfg8 = cfg;
cfg8.method = 'detect';
cfg8.artfctdef.tms.channel = 'MZC01';
cfg8.artfctdef.tms.trlpadding = 0; % cannot be >0 for data in memory
[cfg8, artifact8] = ft_artifact_tms(cfg8, data);

cfg9 = cfg;
cfg9.artfctdef.zvalue.cutoff = 3;
cfg9.artfctdef.zvalue.channel = 'MZC01';
[cfg9, artifact9] = ft_artifact_zvalue(cfg9, data);

%%

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
cfg.continuous = 'yes';
cfg.trl = trlfile;

ft_databrowser(cfg);


