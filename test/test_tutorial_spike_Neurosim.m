function test_tutorial_spike_Neurosim

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_tutorial_spike_sim
% performs all the operations mentioned int the spike tutorial
% (http://www.fieldtriptoolbox.org/tutorial/spike), but only plots figures
% that are called by external functions (e.g. ft_spike_plot_isireturn).
% 
% This corresponds to the tutorial around 25 Sept 2012
% 
% In contrast to test_tutorial_spike this test script uses data generated
% using NeuroSim. 
% (It therefore skips the tutorial part covering waveform analysis)
% 
% This script calls many functions and might take a while to complete.
% TEST ft_read_spike ft_spike_select ft_read_event 
% TEST ft_definetrial ft_spike_maketrials ft_read_header ft_checkdata
% TEST ft_spike_isi ft_spike_plot_isireturn ft_spike_psth ft_spikedensity
% TEST ft_spike_plot_raster ft_spike_xcorr

spike2 = ft_read_spike(dccnpath('/home/common/matlab/fieldtrip/data/test/original/neurosim/spikes')); % should be the folder containing 'signals' and 'spikes'; or the spikes file directly


cfg              = [];
cfg.spikechannel = [1, 4]; % select one inhibitory and one excitatory neuron
spike = ft_spike_select(cfg, spike2);



%% waveform-based cleaning; not useful, or available, with neurosim data
% cfg             = [];
% cfg.fsample     = 40000;
% cfg.interpolate = 1; % keep the density of samples as is
% [wave, spikeCleaned] = ft_spike_waveform(cfg,spike);
% 
% % 
% % for k = [1 2]
% %   figure, 
% %   sl = squeeze(wave.dof(k,:,:))>1000; % only keep samples with enough spikes
% %   plot(wave.time(sl), squeeze(wave.avg(k,:,sl)),'k') % factor 10^6 to get microseconds
% %   hold on
% %  
% %   % plot the standard deviation
% %   plot(wave.time(sl), squeeze(wave.avg(k,:,sl))+sqrt(squeeze(wave.var(k,:,sl))),'k--') 
% %   plot(wave.time(sl), squeeze(wave.avg(k,:,sl))-sqrt(squeeze(wave.var(k,:,sl))),'k--')
% %  
% %   axis tight
% %   set(gca,'Box', 'off') 
% %   xlabel('time')
% %   ylabel('normalized voltage')
% % end

%%


cfg          = []; 
cfg.trialdef.triallength = 0.5; %duration in seconds (can also be 1 or Inf)
cfg.datafile = dccnpath('/home/common/matlab/fieldtrip/data/test/original/neurosim/spikes'); %should be the directory containing both 'spikes' and 'signals'
cfg.trialfun = 'ft_trialfun_general';
cfg = ft_definetrial(cfg);


cfg.trlunit='samples'; %ft_trialfun_general gives us samles, not timestamps
cfg.hdr = spike.hdr; %using 'samples' requires the header information.

% the ft_spike functions use pedantic cfg checking, so remove the trialdef
% field which is illegal here
cfg = rmfield(cfg, 'trialdef');
spikeTrials = ft_spike_maketrials(cfg,spike); 


dat = ft_checkdata(spikeTrials,'datatype', 'raw', 'fsample', 1000);


cfg       = [];
cfg.bins  = [0:0.0005:0.1]; % use bins of 0.5 milliseconds
cfg.bins  = [0:0.0001:0.1]; % use bins of 0.1 milliseconds
cfg.param = 'coeffvar'; % compute the coefficient of variation (sd/mn of isis)
isih = ft_spike_isi(cfg,spikeTrials);

for k = [1 2] % only do for the single units
  cfg              = [];
  cfg.spikechannel = isih.label{k};
  cfg.interpolate  = 5; % interpolate at 5 times the original density
  cfg.window       = 'gausswin'; % use a gaussian window to smooth
  cfg.winlen       = 0.004; % the window by which we smooth has size 4 by 4 ms
  cfg.colormap     = jet(300); % colormap
  cfg.scatter      = 'no'; % do not plot the individual isis per spike as scatters
  figure, ft_spike_plot_isireturn(cfg,isih);
end


 


% cfg             = []; 
% cfg.binsize     =  0.1; % if cfgPsth.binsize = 'scott' or 'sqrt', we estimate the optimal bin size from the data itself
% cfg.outputunit  = 'rate'; % give as an output the firing rate
% cfg.latency     = [-1 3]; % between -1 and 3 sec.
% cfg.vartriallen = 'yes'; % variable trial lengths are accepted
% cfg.keeptrials  = 'yes'; % keep the psth per trial in the output
% psth = ft_spike_psth(cfg,spikeTrials);

cfg         = [];
cfg.binsize =  [0.05];
cfg.latency = [0 0.4];
psth = ft_spike_psth(cfg,spikeTrials);
 
cfg              = [];
cfg.topplotfunc  = 'line'; % plot as a line
cfg.spikechannel = spikeTrials.label([1 2]);
cfg.latency      = [-1 3];
cfg.errorbars    = 'std'; % plot with the standard deviation
cfg.interactive  = 'no'; % toggle off interactive mode
figure, ft_spike_plot_raster(cfg,spikeTrials, psth); 


cfg         = [];
cfg.latency = [0 0.4]; 
cfg.timwin  = [-0.025 0.025];
cfg.fsample = 1000; % sample at 1000 hz
sdf = ft_spikedensity(cfg,spikeTrials);
cfg              = [];
cfg.topplotfunc  = 'line'; % plot as a line plot
cfg.spikechannel = spikeTrials.label([1 2]); % can also select one unit here
cfg.latency      = [-1 3];
cfg.errorbars    = 'std'; % plot with standard deviation
cfg.interactive  = 'no'; % toggle off interactive mode
figure, ft_spike_plot_raster(cfg,spikeTrials, sdf) 


 
cfg              = [];
cfg.spikechannel = [1, 2, 4, 5]; % select two inhibitory and two excitatory neurons
spike = ft_spike_select(cfg, spike2);
 
cfg          = []; 
cfg.trialdef.triallength = 0.5; %duration in seconds (can also be 1 or Inf)
cfg.datafile = dccnpath('/home/common/matlab/fieldtrip/data/test/original/neurosim/spikes'); %should be the directory containing both 'spikes' and 'signals'
cfg.trialfun = 'ft_trialfun_general';
cfg = ft_definetrial(cfg);


cfg.trlunit='samples'; %ft_trialfun_general gives us samles, not timestamps
cfg.hdr = ft_read_header(cfg.datafile);
% the ft_spike functions use pedantic cfg checking, so remove the trialdef
% field which is illegal here
cfg = rmfield(cfg, 'trialdef');
spikeTrials = ft_spike_maketrials(cfg,spike); 

cfg             = [];
cfg.maxlag      = 0.2; % maximum 200 ms
cfg.binsize     = 0.001; % bins of 1 ms
cfg.outputunit  = 'proportion'; % make unit area
cfg.latency     = [0 0.3];
cfg.vartriallen = 'no'; % do not allow variable trial lengths
cfg.method      = 'xcorr'; % compute the normal cross-correlogram
Xc = ft_spike_xcorr(cfg,spikeTrials);  
 
% compute the shuffled correlogram
cfg.method      = 'shiftpredictor'; % compute the shift predictor
Xshuff = ft_spike_xcorr(cfg,spikeTrials);

% iCmb = 3;
% jCmb = 4;
% figure
% xcSmoothed = conv(squeeze(Xc.xcorr(iCmb,jCmb,:)),ones(1,5)./5,'same'); % do some smoothing
% hd = plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'k'); % leave out borders (because of smoothing)
% hold on
% xcSmoothed = conv(squeeze(Xshuff.shiftpredictor(iCmb,jCmb,:)),ones(1,5)./5,'same');    
% plot(Xc.time(3:end-2),xcSmoothed(3:end-2),'r')
% hold on
% xlabel('delay')
% ylabel('proportion of coincidences')        
% title([Xc.label{iCmb} Xc.label{jCmb}])
% axis tight

% compute the spike densities
cfg         = [];
cfg.latency = [0 0.3]; 
cfg.timwin  = [-0.035 0.025];
cfg.fsample = 200;
[sdf] = ft_spikedensity(cfg,spikeTrials);
 
% compute the joint psth
cfg               = [];
cfg.normalization = 'no';
cfg.channelcmb    = spikeTrials.label(3:4);
cfg.method        = 'jpsth';
jpsth = ft_spike_jpsth(cfg,sdf);
 
cfg.method = 'shiftpredictor';
jpsthShuff = ft_spike_jpsth(cfg,sdf);
%  
% % subtract the predictor
% jpsthSubtr = jpsth;
% jpsthSubtr.jpsth = jpsth.jpsth-jpsthShuff.shiftpredictor;
% 
% 
% cfg        = [];
% figure
% ft_spike_plot_jpsth(cfg,jpsth)
% figure
% ft_spike_plot_jpsth(cfg,jpsthShuff)
% figure
% ft_spike_plot_jpsth([],jpsthSubtr)
end
