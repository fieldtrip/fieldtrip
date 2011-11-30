function test_ft_networkanalysis

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

cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
stat = ft_networkanalysis(cfg, tmp);