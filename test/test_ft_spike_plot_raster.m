function test_ft_spike_plot_isireturn()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spike_plot_raster
% TEST test_ft_spike_plot_raster

%%
clear
spikesPerTrial = 10;
nTrials = 40;
shapePar = 2;
scalePar = 3;
time       = linspace(0,1,1000);
data.trial(1:nTrials) = {zeros(1,length(time))};
data.time(1:nTrials) = {time};
for iTrial = 1:nTrials    
  for iUnit = 1:3
      if iUnit==1
        data.trial{iTrial}(iUnit,:) = rand(1,length(time));
      else
        spikeTimes = 0.015*gamrnd(shapePar*iUnit,scalePar,[spikesPerTrial 1]) + 0.001*iTrial;
        spikeTimes(spikeTimes>1) = [];  
        smp = [];
        for iSpike = 1:length(spikeTimes)
          smp(iSpike)        = nearest(time,spikeTimes(iSpike));
        end
        data.trial{iTrial}(iUnit,smp) = 1;
      end
  end
end
data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
data.label{1} = 'chan1';
data.label{end+1} = 'spk2';
data.label{end+1} = 'spk3';

%% create the spike format as well by means of our conversion script
% show that the psth works also with the poisson format
spike = ft_checkdata(data,'datatype', 'spike', 'feedback', 'yes');

%% create the first rasterplot, using both types of data format
% here, we auto-color the different units

cfgRaster = [];
cfgRaster.spikelength = 0.8;
cfgRaster.linewidth = 1;
cfgRaster.spikechannel = [];
figure
H = ft_spike_plot_raster(cfgRaster,spike);
%%

cfgRaster = [];
cfgRaster.spikelength = 0.8;
cfgRaster.linewidth = 1;
cfgRaster.trials = 2:2:20;
figure
H = ft_spike_plot_raster(cfgRaster,spike);


%% we compute the psth by calling the psth function
psthData = ft_spike_psth([],spike);
%% demonstrate plotting of psth with our rasterplot
figure
cfgRaster =[];
cfgRaster.spikechannel = [1 2];
cfgRaster.topplotsize = 0.2;
H = ft_spike_plot_raster(cfgRaster,spike, psthData);





