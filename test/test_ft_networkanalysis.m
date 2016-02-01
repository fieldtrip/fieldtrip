function test_ft_networkanalysis

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_networkanalysis ft_networkanalysis

data = [];
for k = 1:5
  x = randn(7,15);
  cx = x*x';
  cx = cx./sqrt(diag(x)*diag(x)');
  data.cohspctrm(:,:,k) = cx;
end
data.freq  = 1:5;
data.label = {'chan1';'chan2';'chan3';'chan4';'chan5';'chan6';'chan7'};
data.dimord = 'chan_chan_freq';
data.cfg    = 'this is the cfg';

tmp = data;
tmp.cohspctrm = data.cohspctrm>0.3;

% at present just checks for undirected binary and weighted graphs
cfg           = [];
cfg.parameter = 'cohspctrm';

cfg.method    = 'assortativity';
stat1 = ft_networkanalysis(cfg, tmp);
stat2 = ft_networkanalysis(cfg, data);

cfg.method    = 'betweenness';
stat3 = ft_networkanalysis(cfg, tmp);
stat4 = ft_networkanalysis(cfg, data);

cfg.method    = 'clustering_coef';
stat5 = ft_networkanalysis(cfg, tmp);
stat6 = ft_networkanalysis(cfg, data);

cfg.method    = 'degrees';
stat7 = ft_networkanalysis(cfg, tmp);
stat8 = ft_networkanalysis(cfg, data);
