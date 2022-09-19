function test_issue1866

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_checkdata ft_connectivityplot


data.trial{1} = randn(2,1000);
data.trial{2} = randn(2,1000);
data.trial{3} = randn(2,1000);
data.time{1}=(0:999)./1000-0.5;
data.time{2}=(0:999)./1000-0.5;
data.time{3}=(0:999)./1000-0.5;
data.label={'a';'b'};


cfg = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:1:30;
cfg.t_ftimwin    = 4./cfg.foi;  % 4 cycles per time window
cfg.toi          = -0.5:0.05:0.6;
cfg.channel      = 'all';

f = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'coh';
c = ft_connectivityanalysis(cfg, f);


cfg = [];
cfg.parameter = 'cohspctrm';
figure; ft_connectivityplot(cfg, c);

cfg.latency = 0;

% this does not 'work', i.e. ft_selectdata is not executed
figure; ft_connectivityplot(cfg, c);

% this should work
figure; ft_connectivityplot(cfg, ft_selectdata(cfg, c));

f2 = ft_checkdata(f, 'cmbstyle', 'full');
assert(isequal(f.freq, f2.freq));

