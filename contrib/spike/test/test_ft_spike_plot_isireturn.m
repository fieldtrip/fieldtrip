function test_ft_spike_plot_isireturn()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spike_plot_isireturn
% TEST ft_spike_plot_isireturn

nTrials = 100;
data = [];
time       = linspace(0,1,1000);
data.trial(1:nTrials) = {zeros(1,length(time))};
data.time(1:nTrials) = {time};
for iTrial = 1:nTrials    
  for iUnit = 1:3
      data.trial{iTrial}(iUnit,:) = double(rand(1,1000)<0.05);
  end
end
data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
data.label{1} = 'chan1';
data.label{end+1} = 'chan2';
data.label{end+1} = 'chan3';

% show that the psth works also with the poisson format
spike = ft_checkdata(data,'datatype', 'spike', 'feedback', 'yes');
for iUnit  = 1:3
  spike.time{iUnit} =   spike.time{iUnit} + 0.001*rand(1,length(spike.time{iUnit}))
end
%%
% now test the isi
cfgIsi = [];
cfgIsi.keeptrials = 'yes';
cfgIsi.param='gamfit';
isih = ft_spike_isi(cfgIsi,spike);
%%
figure
cfgIsi = [];
cfgIsi.spikechannel = 1;
cfgIsi.plotfit = 'yes'
cfgIsi.latency = [0 0.2];
h = ft_spike_plot_isi(cfgIsi,isih)
%%
cfgIsi = [];
cfgIsi.keeptrials = 'yes';
isih = ft_spike_isi(cfgIsi,spike);
%%
cfgIsi = [];
cfgIsi.spikechannel = 1;
isihS = ft_spike_isi(cfgIsi,spike);
%%
cfgIsi = [];
cfgIsi.spikechannel = 'all';
isihS = ft_spike_isi(cfgIsi,spike);
%%
cfgRet = [];
cfgRet.window = 'gausswin';
cfgRet.winlen = 0.2
cfgRet.interpolate = 1;
cfgRet.scattersize = 1;
cfgRet.density = 'yes'
cfgRet.spikechannel = 1;
figure
H = ft_spike_plot_isireturn(cfgRet,isihS);
%%
figure
H = ft_spike_plot_isireturn([],isihS);
%%

