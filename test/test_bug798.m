function test_bug798

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqstatistics ft_selectdata ft_datatype_freq ft_appendfreq

% note that this bug is related to bug 921

datadir = (dccnpath('/home/common/matlab/fieldtrip/data/test/bug798'));
load(fullfile(datadir,'t2_subj1.mat'));
load(fullfile(datadir,'t2_subj1_null'));
load(fullfile(datadir,'t2_subj2'));
load(fullfile(datadir,'t2_subj2_null'));

% the following gave an error on 6 July 2011, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=798
a = ft_checkdata(t2_subj1, 'datatype', 'freq');

% after changing the ft_datatype and ft_checkdata the following also works
a = ft_checkdata(t2_subj1, 'datatype', 'timelock');
% in both cases it adds a dimord and a fake (nan) dimension for time/freq

% the reported problem was something like this, resulting in 274x0 sized outputs
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_pooledT';
cfg.tail = 0;
cfg.alpha = 0.025;
cfg.correctm = 'cluster';
cfg1 = [];
cfg1.gradfile = fullfile(datadir,'ctf275.mat');
cfg1.method = 'triangulation';
cfg1.feedback = 'yes';
cfg.neighbours = ft_prepare_neighbours(cfg1, t2_subj1);
cfg.numrandomization = 10;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric';
cfg.clusteralpha     = 0.05;
cfg.clustercritval   = [-1.96 1.96];
cfg.clustertail      =  0;
cfg.design(1,1:2*2) = [ones(1,2) 2*ones(1,2)];
cfg.design(2,1:2*2) = [1:2 1:2];
cfg.ivar =1;
cfg.uvar =2;
stat = ft_freqstatistics(cfg,t2_subj1,t2_subj2,t2_subj1_null,t2_subj2_null);
assert(all(size(stat.prob)==[274 1]));

% the bug reduces to
data = ft_checkdata(t2_subj1, 'datatype', 'freq');

tmpcfg = [];
tmpcfg.avgoverchan = 'no';
tmpcfg.frequency = [-inf inf];
tmpcfg.avgoverfreq = 'no';
data = ft_selectdata(tmpcfg, data);
% which is due to ft_checkdata inserting a nan in the fake frequency axis
% fixed by roboos on 1 September 2011
assert(all(size(data.powspctrm)==[274 1]));

% then the bug is in
argin{1} = ft_checkdata(t2_subj1, 'datatype', 'freq');
argin{2} = ft_checkdata(t2_subj1_null, 'datatype', 'freq');
argin{3} = ft_checkdata(t2_subj2, 'datatype', 'freq');
argin{4} = ft_checkdata(t2_subj2_null, 'datatype', 'freq');
% and

% but this should now work
cfg = [];
cfg.parameter = 'powspctrm';
cfg.appenddim = 'rpt';
data2 = ft_appendfreq(cfg, argin{:});
assert(~isfield(data2, 'stat'));
assert(~isfield(data2, 'prob'));
assert(strcmp(data2.dimord, 'rpt_chan_freq'));
