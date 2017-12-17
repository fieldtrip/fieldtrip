function test_bug1483

% MEM 1500mb
% WALLTIME 00:10:00

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
dataold1 = ft_selectdata(data, 'rpt', 2:4);
dataold2 = ft_selectdata(data, 'channel', data.label(1));

freqold1 = ft_selectdata(freq, 'rpt', 2:4);
freqold2 = ft_selectdata(freq, 'channel', freq.label(1));


% new style
cfg = [];
cfg.trials = 2:4;
datanew1 = ft_selectdata(cfg, data);
freqnew1 = ft_selectdata(cfg, freq);

cfg = [];
cfg.channel = data.label(1);
datanew2 = ft_selectdata(cfg, data);
freqnew2 = ft_selectdata(cfg, freq);


% BUG CONFIRMED

%%



