function test_bug1483

% TEST test_bug1483
% TEST ft_selectdata ft_selectdata_new ft_selectdata_old

%% first confirm the reported bug, i.e. the incapability of ft_selectdata_new
% to select 'rpt'

% create some dummy data
for k = 1:5
  data.trial{k} = ones(2,3).*k;
  data.time{k}  = -1:1;
end
data.label = {'chan1';'chan2'};

freq.powspctrm = rand(5,2,3);
freq.dimord    = 'rpt_chan_freq';
freq.label     = {'chan1';'chan2'};
freq.freq      = [1 2 3];
freq.cumtapcnt = ones(5,1);

% old style
dataold = ft_selectdata(data, 'rpt', 2:4);
freqold = ft_selectdata(freq, 'rpt', 2:4);

% new style
cfg = [];
cfg.trials = 2:4;
datanew = ft_selectdata(cfg, data);
freqnew = ft_selectdata(cfg, freq);

% BUG CONFIRMED

%%



