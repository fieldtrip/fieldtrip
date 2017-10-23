function test_tutorial_connectivity20130308

% MEM 1500mb
% WALLTIME 00:10:00


% Simulated data with directed connections
% We will first simulate some data with a known connectivity structure built in. This way we know what to expect in terms of connectivity. To simulate data we use ft_connectivitysimulation. We will use an order 2 multivariate autoregressive model. The necessary ingredients are a set of NxN coefficient matrices, one matrix for each time lag. These coefficients need to be stored in the cfg.param field. Next to the coefficients we have to specify the NxN covariance matrix of the innovation noise. This matrix needs to be stored in the cfg.noisecov field. The model we are going to use to simulate the data is as follows:
% 
% x(t) = 0.8*x(t-1) - 0.5*x(t-2)
% 
% y(t) = 0.9*y(t-1) + 0.5*z(t-1) - 0.8*y(t-2)
% 
% z(t) = 0.5*z(t-1) + 0.4*x(t-1) - 0.2*z(t-2)

cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';

cfg.params(:,:,1) = [ 0.8    0    0 ; 
                        0  0.9  0.5 ;
                      0.4    0  0.5];
                      
cfg.params(:,:,2) = [-0.5    0    0 ; 
                        0 -0.8    0 ; 
                        0    0 -0.2];
                        
cfg.noisecov      = [ 0.3    0    0 ;
                        0    1    0 ;
                        0    0  0.2];

data              = ft_connectivitysimulation(cfg);

% The simulated data consists of 3 channels in 500 trials. You can easily visualize the data for example in the first trial using
% 
figure
plot(data.time{1}, data.trial{1}) 
legend(data.label)
xlabel('time (s)')


% or browse through the complete data using
% 
cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly' 
ft_databrowser(cfg, data);

% 
% 
% Computation of the multivariate autoregressive model
% To be able to compute spectrally resolved Granger causality, or other frequency-domain directional measures of connectivity, we have to fit an autoregressive model to the data. This is done using the ft_mvaranalysis function.
% 
% For the actual computation of the autoregressive coefficients FieldTrip makes use of an implementation from third party toolboxes. At present ft_mvaranalysis supports the biosig and bsmart toolboxes for these computations.
% 
% In this tutorial we will use the bsmart toolbox. The relevant functions have been included in the FieldTrip release in the fieldtrip/external/bsmart directory.
% 
cfg         = [];
cfg.order   = 5;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data);

% mdata = 
%          dimord: 'chan_chan_lag'
%           label: {3x1 cell}
%          coeffs: [3x3x5 double]
%        noisecov: [3x3 double]
%             dof: 500
%     fsampleorig: 200
%             cfg: [1x1 struct]
%             
% The resulting variable mdata contains a description of the data in terms of a multivariate autoregressive model. For each time-lag up to the model order (which is 5 in this case), a 33 matrix of coefficients is outputted. The noisecov-field contains covariance matrix of the model's residuals.
% 
% Exercise 1
% Compare the parameters specified for the simulation with the estimated coefficients and discuss.
% 
% Computation of the spectral transfer function
% From the autoregressive coefficients it is now possible to compute the spectral transfer matrix, for which we use ft_freqanalysis.
% 
cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);
 
% mfreq = 
%         label: {3x1 cell}
%          freq: [1x101 double]
%        dimord: 'chan_chan_freq'
%      transfer: [3x3x101 double]
%      noisecov: [3x3 double]
%     crsspctrm: [3x3x101 double]
%           dof: 500
%           cfg: [1x1 struct]
%  
% The resulting mfreq data structure contains the pairwise transfer function between the 3 channels for 101 frequencies.
% 
% It is also possible to compute the spectral transfer function using non-parametric spectral factorization of the cross-spectral density matrix. For this, we need a Fourier decomposition of the data. This is done in the following section.
% 
% 
% Non-parametric computation of the cross-spectral density matrix
% Some connectivity metrics can be computed from a non-parametric spectral estimate (i.e. after the application of the FFT-algorithm and conjugate multiplication to get cross-spectral densities), such as coherence, phase-locking value and phase slope index. The following part computes the fourier-representation of the data using ft_freqanalysis. It is not necessary to compute the cross-spectral density at this stage, because the function used in the next step, ft_connectivityanalysis, contains functionality to compute the cross-spectral density from the fourier coefficients.
% 
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, data);

% freq = 
%             label: {3x1 cell}
%            dimord: 'rpttap_chan_freq'
%              freq: [1x101 double]
%     fourierspctrm: [1500x3x101 double]
%         cumsumcnt: [500x1 double]
%         cumtapcnt: [500x1 double]
%               cfg: [1x1 struct]
% The resulting freq structure contains the spectral estimate for 3 tapers in each of the 500 trials (hence 1500 estimates), for each of the 3 channels and for 101 frequencies.
% 
% 
% Computation and inspection of the connectivity measures
% The actual computation of the connectivity metric is done by ft_connectivityanalysis. This function is transparent to the type of input data, i.e. provided the input data allows the requested metric to be computed, the metric will be calculated. Here, we provide an example for the computation and visualization of the coherence coefficient.
% 
cfg           = [];
cfg.method    = 'coh';
coh           = ft_connectivityanalysis(cfg, freq);
cohm          = ft_connectivityanalysis(cfg, mfreq);
% Subsequently, the data can be visualized using ft_connectivityplot.
% 
cfg           = [];
cfg.parameter = 'cohspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, coh, cohm);

% 
% The coherence measure is a symmetric measure, which means that it does not provide information regarding the direction of information flow between any pair of signals. In order to analyze directionality in interactions, measures based on the concept of granger causality can be computed. These measures are based on an estimate of the spectral transfer matrix, which can be computed in a straightforward way from the multivariate autoregressive model fitted to the data.
% 
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, granger);


% Instead of plotting it with ft_connectivityplot, you can use the following low-level MATLAB plotting code which gives a better understanding of the numerical representation of the results.
% 
figure
for row=1:3
for col=1:3
  subplot(3,3,(row-1)*3+col);
  plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
  ylim([0 1])
end
end

% 
% Exercise 2
% Discuss the differences between the granger causality spectra, and the coherence spectra.
% Exercise 3
% Compute the following connectivity measures from the mfreq data, and visualize and discuss the results: partial directed coherence (pdc), directed transfer function (dtf), phase slope index (psi)
% 
% Simulated data with common pick-up and different noise levels
%  this is under progress
% 
% When working with electrophysiological data (EEG/MEG/LFP) the signals that are picked up by the individual channels invariably consist of instantaneous mixtures of the underlying source signals. This mixing can severely affect the outcome of connectivity analysis, and thus affects the interpretation. We will demonstrate this by simulating data in 2 channels, where each of the channels consists of a weighted combination of temporally white noise unique to each of the channels, and a common input of a band-limited signal (filtered between 15 and 25 Hz). We will compute connectivity between these channels, and show that the common input can give rise to spurious estimates of connectivity.
% 
% % create some instantaneously mixed data
% 
% define some variables locally
nTrials  = 100;
nSamples = 1000;
fsample  = 1000;

% mixing matrix
mixing   = [0.8 0.2 0;
              0 0.2 0.8];

data       = [];
data.trial = cell(1,nTrials);
data.time  = cell(1,nTrials);
for k = 1:nTrials
  dat = randn(3, nSamples);
  dat(2,:) = ft_preproc_bandpassfilter(dat(2,:), 1000, [15 25]);
  dat = 0.2.*(dat-repmat(mean(dat,2),[1 nSamples]))./repmat(std(dat,[],2),[1 nSamples]);
  data.trial{k} = mixing * dat;
  data.time{k}  = (0:nSamples-1)./fsample;
end
data.label = {'chan1' 'chan2'}';

figure;plot(dat'+repmat([0 1 2],[nSamples 1]));
title('original ''sources''');

figure;plot((mixing*dat)'+repmat([0 1],[nSamples 1])); 
axis([0 1000 -1 2]);
set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
title('mixed ''sources''');
 

% do spectral analysis
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'fourier';
cfg.foilim    = [0 200];
cfg.tapsmofrq = 5;
freq          = ft_freqanalysis(cfg, data);
fd            = ft_freqdescriptives(cfg, freq);

figure;plot(fd.freq, fd.powspctrm);
set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
title('powerpectrum');

% 
% compute connectivity
cfg = [];
cfg.method = 'granger';
g = ft_connectivityanalysis(cfg, freq);
cfg.method = 'coh';
c = ft_connectivityanalysis(cfg, freq);
% visualize the results
cfg = [];
cfg.parameter = 'grangerspctrm';
ft_connectivityplot(cfg, g);
cfg.parameter = 'cohspctrm';
ft_connectivityplot(cfg, c);


% Exercise 4
% Simulate new data using the following mixing matrix:
% [0.9 0.1 0;0 0.2 0.8] 
% and recompute the connectivity measures. Discuss what you see.
% 
% Exercise 5
% Play a bit with the parameters in the mixing matrix and see what is the effect on the estimated connectivity.
% Exercise 6
% Simulate new data where the 2 mixed signals are created from 4 underlying sources, and where two of these sources are common input to both signals, and where these two sources are temporally shifted copies of one another.
% Hint: the mixing matrix could look like this:
% 
% [a b c 0; 0 d e f];
% and the trials could be created like this:
% 
% for k = 1:nTrials
%   dat = randn(4, nSamples+10);
%   dat(2,:) = ft_preproc_bandpassfilter(dat(2,:), 1000, [15 25]);
%   dat(3,1:(nSamples)) = dat(2,11:(nSamples+10)); 
%   dat = dat(:,1:1000);
%   dat = 0.2.*(dat-repmat(mean(dat,2),[1 nSamples]))./repmat(std(dat,[],2),[1 nSamples]);
%   data.trial{k} = mixing * dat;
%   data.time{k}  = (0:nSamples-1)./fsample;
% end
% Compute connectivity between the signals and discuss what you observe. In particular, also compute measures of directed interaction.
% 
% 
% Connectivity between MEG virtual channel and EMG
% The previous two examples were using simulated data, either with a clear directed connectivity structure, or with a trivial pick-up of a common source in two channels. We will now continue with connectivity analysis on real MEG data. The dataset is the same as the one used in the Analysis of corticomuscular coherence tutorial.
% 
% In short, the dataset consists of combined MEG and EMG recordings while the subject lifted his right hand. The coherence tutorial introduction contains a more elaborate description of the experiment and the dataset and a detailed analysis can be found in the corresponding paper 1). Due to the long distance between the EMG and the MEG, there is no volume conduction and hence no common pick-up. Hence this dataset lends itself well for connectivity analysis. But rather than doing an analysis between the EMG and one of the MEG channels (as in the original study), we will extract the cortical activity using a beamformer virtual channel.
% 
% 
% Compute the spatial filter for the region of interest
% We start with determining the motor cortex as region of interest. At the end of the coherence tutorial it is demonstrated how to make a 3-D reconstruction of the cortico-muscolar coherence (CMC) using the DICS algorithm. That source reconstruction serves as starting point for this analysis.
% 
% You can download the result from the DICS reconstruction from the FieldTrip ftp server (source.mat)
% 
% We will first determine the position on which the cortico-muscular coherence is the largest.
% 

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/connectivity'));
load source

[maxval, maxindx] = max(source.avg.coh);
maxpos = source.pos(maxindx,:);

% maxpos = 
%     4 -3 12
% The cortical position is expressed in individual subject head-coordinates and in centimeter. Relative to the center of the head (in between the ears) the position is 4 cm towards the nose, -3 towards the left side (i.e., 3 cm towards the right!) and 12 cm towards the vertex.
% 
% The ft_sourceanalysis methods are usually applied to the whole brain using a regular 3-D grid or using a triangulated cortical sheet. You can also just specify the location of a single or multiple points of interest with cfg.grid.pos and the LCMV beamformer will simply be performed at the location of interest.
% 
% The LCMV beamformer spatial filter for the location of interest will pass the activity at that location with unit-gain, while optimally suppressing all other noise and other source contributions to the MEG data. The LCMV implementation in FieldTrip requires the data covariance matrix to be computed with ft_timelockanalysis.
% 
% Rather than doing all the preprocessing again, you can download the preprocessed data from the FieldTrip ftp server (data.mat)
% 
load data
% 
%% compute the beamformer filter
cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
timelock              = ft_timelockanalysis(cfg, data);

cfg             = [];
cfg.method      = 'lcmv';
cfg.hdmfile     = 'SubjectCMC.hdm';
cfg.grid.pos    = maxpos;
cfg.keepfilter  = 'yes';
source          = ft_sourceanalysis(cfg, timelock);
% The source reconstruction contains the estimated power and the source-level time-series of the averaged ERF, but here we are not interested in those. The cfg.keepfilter option results in the spatial filter being kept in the output source structure. That spatial can be used to reconstruct the single-trial time series as a virtual channel by multiplying it with the original MEG data.
% 
% 
% Extract the virtual channel time-series
%% construct the 3-D virtual channel at the location of interest
beamformer = source.avg.filter{1};

chansel = ft_channelselection('MEG', data.label); % find the names
chansel = match_str(data.label, chansel);         % find the indices

sourcedata = [];
sourcedata.label = {'x', 'y', 'z'};
sourcedata.time = data.time;
for i=1:length(data.trial)
  sourcedata.trial{i} = beamformer * data.trial{i}(chansel,:);
end
% The LCMV spatial filter is computed here without applying any time-domain filters. Consequently, it will have to suppress all noise in the data in all frequency bands. The spatial filter derived from the broadband data allows us to compute a broadband source level time-series.
% If you would know that the subsequent analysis would be limited to a specific frequency range in the data (e.g. everything above 30 Hz), you could first apply a filter using ft_preprocessing (e.g. cfg.hpfilter=yes and cfg.hpfreq=30) prior to computing the covariance and the spatial filter.
% 
% The sourcedata structure resembles the raw-data output of ft_preprocessing and consequently can be used in any follow-up function. You can for example visualize the single-trial virtual channel time-series using ft_databrowser:
% 
cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
ft_databrowser(cfg, sourcedata);

% 
% Notice that the reconstruction contains three channels, for the x-, the y- and the z-component of the equivalent current dipole source at the location of interest.
% 
% 
% Project along the strongest dipole direction
% The interpretation of connectivity is facilitated if we can compute it between two plain channels rather than between one channel versus a triplet of channels. Therefore we will project the time-series along the dipole direction that explains most variance. This projection is equivalent to determining the largest (temporal) eigenvector and can be computationally performed using the singular value decomposition (svd).
% 
%% construct a single virtual channel in the maximum power orientation
timeseries = cat(2, sourcedata.trial{:});

[u, s, v] = svd(timeseries, 'econ');

% whos u s v
%   Name           Size              Bytes  Class     Attributes
% 
%   s              3x3                  72  double              
%   u              3x3                  72  double              
%   v         196800x3             4723200  double            
% Matrix u contains the spatial decomposition, matrix v the temporal and on the diagonal of matrix s you can find the eigenvalues. See "help svd" for more details.
% 
% We now recompute the virtual channel time-series, but now only for the dipole direction that has the most power.
% 
% this is equal to the first column of matrix V, apart from the scaling with s(1,1)
timeseriesmaxproj = u(:,1)' * timeseries;
virtualchanneldata = [];
virtualchanneldata.label = {'cortex'};
virtualchanneldata.time = data.time;
for i=1:length(data.trial)
  virtualchanneldata.trial{i} = u(:,1)' * beamformer * data.trial{i}(chansel,:);
end
% Rather than using a sourcemodel in the beamformer that consists of all three (x, y, z) directions, you can also have the beamformer compute the filter for only the optimal source orientation. This is implemented using the cfg.lcmv.fixedori='yes' option.
% Recompute the spatial filter for the optimal source orientation and using that spatial filter (a 1151 vector) recompute the time-series.
% 
% Investigate and describe the difference between the two time-series. What is the difference between the two dipole orientations?
% 
% Note that one orientation is represented in the SVD matrix "u" and the other is in the source.avg.ori field.
% 
% 
% Combine the virtual channel with the EMG
% The raw data structure containing one (virtual) channel can be combined with the two EMG channels from the original preprocessed data.
% 
%% select the two EMG channels
cfg = [];
cfg.channel = 'EMG';
emgdata = ft_selectdata(cfg, data);

%% combine the virtual channel with the two EMG channels
cfg = [];
combineddata = ft_appenddata(cfg, virtualchanneldata, emgdata);

% save combineddata combineddata
% 
% Compute the connectivity
% The resulting combined data structure has three channels: the activity from the cortex, the left EMG and the right EMG. We can now continue with regular channel-level connectivity analysis.
% 
%% compute the spectral decomposition
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'cortex' 'EMGlft' 'EMGrgt'};
freq    = ft_freqanalysis(cfg, combineddata);

cfg = [];
cfg.method = 'coh';
coherence = ft_connectivityanalysis(cfg, freq);
% This computes the spectral decomposition and the coherence spectrum between all channel pairs, which can be plotted with
% 
cfg = [];
cfg.zlim = [0 0.2];
figure
ft_connectivityplot(cfg, coherence);
title('coherence')

% 
% To look in more detail into the numerical representation of the coherence results, you can use
% 
figure
plot(coherence.freq, squeeze(coherence.cohspctrm(1,2,:)))
title(sprintf('connectivity between %s and %s', coherence.label{1}, coherence.label{2}));
xlabel('freq (Hz)')
ylabel('coherence')

% 
% The spectrum reveals coherence peaks at 10 and 20 Hz (remember that the initial DICS localizer was done at beta). Furthermore, there is a broader plateau of coherence in the gamma range from 40-50 Hz.
% 
% The spectral decomposition was performed with mutitapering and 5 Hz spectral smoothing (i.e. 5Hz in both directions). Recompute the spectral decomposition and the coherence with a hanning taper. Recompute it with mutitapering and 10 Hz smoothing. Plot the three coherence spectra and look at the differences.
% Rather than looking at undirected coherence, the virtual channel level data can now also easily be submitted to directed connectivity measures. Compute the spectrally resolved granger connectivity and try to assess whether the directionality is from cortex to EMG or vice versa.
% 
% Summary and further reading
% This tutorial demonstrates how to compute connectivity measures between two time series. If you want to learn how to make a distributed representation of connectivity throughout the whole brain, you may want to continue with the corticomuscular coherence tutorial.
% 
% This tutorial was last tested by Robert with revision 6026 of FieldTrip (~20120611) on a 64-bit Mac OS X machine using MATLAB 2011b.
