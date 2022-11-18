function test_tutorial_networkanalysis_eeg20220126

% WALLTIME 00:30:00
% MEM 10gb
% DEPENDENCY ft_networkanalysis

%% # Preprocessing
%
% ## Reading the data

% load EEG data
datadir = '/home/common/matlab/fieldtrip/data/ftp/tutorial/networkanalysis_eeg';
load(fullfile(datadir, 'data_eeg_reref_ica.mat'));
load(fullfile(datadir, 'elec.mat'));

% select EEG electrodes only
cfg         = [];
cfg.channel = elec.label;
data        = ft_selectdata(cfg,data_eeg_reref_ica);
data        = rmfield(data,'grad');
% convert elec positions in mm
elec        = ft_convert_units(elec,'mm');
data.elec   = elec;

%% ## Prepare electrode layout for plotting
% prepare layout and plot
cfg         = [];
cfg.elec    = elec;
layout      = ft_prepare_layout(cfg);
% scale the layout to fit the head outline
lay         =layout;
lay.pos     =layout.pos./.7;
lay.pos(:,1)=layout.pos(:,1)./.9;
lay.pos(:,2)=layout.pos(:,2)+.08;
lay.pos(:,2)=lay.pos(:,2)./.7;
figure;
ft_plot_layout(lay)

%% ## Data segmentation
%
% Next, the data is segmented into overlapping segemnts of 1 second length.
% resegment the data into 1 sec chunks
cfg         = [];
cfg.length  = 1;
cfg.overlap = .5;
dataseg     = ft_redefinetrial(cfg,data);

%% # Spectral analysis and peak picking
% compute the power spectrum
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.keeptrials   = 'no';
datapow          = ft_freqanalysis(cfg, dataseg);
% plot the topography and the spectrum
figure(1);

cfg             = [];
cfg.layout      = lay;
cfg.xlim        = [9 11];
subplot(1,2,1); ft_topoplotER(cfg, datapow);

cfg             = [];
cfg.channel     = {'EEG087', 'EEG088'};
cfg.xlim        = [3 30];
subplot(1,2,2); ft_singleplotER(cfg, datapow);

%% # Computation of the forward model
% load the required geometrical information
load(fullfile(datadir, 'dkatlas.mat'));
load(fullfile(datadir, 'headmodel_eeg.mat'));
load(fullfile(datadir, 'sourcemodel.mat'));

%% visualize the coregistration of sensors, headmodel, and sourcemodel.
figure(2);
% make the headmodel surface transparent
ft_plot_headmodel(headmodel_eeg, 'edgecolor', 'none'); alpha 0.4
ft_plot_sens(dataseg.elec);
view([45 -15 0])

% In Figure 3 it is apparent that the electrodes do not align with the scalp surface. 
% To achieve this we use ft_electroderealign in an interactive mode. 
% Figure 4 provides the settings that had been used to align the electrodes. 
% In particular, the option rotate, scale and translate in Figure 3.
%%
% cfg         = [];
% cfg.method  = 'interactive';
% cfg.headshape = headmodel_eeg.bnd(1);
% cfg.elec    = elec;
% elec_aligned = ft_electroderealign(cfg);

load(fullfile(datadir, 'elec_aligned'));

% make sure the aligned electrodes are updated
dataseg.elec = elec_aligned;

% Before we proceed it is always useful to check the corregistration between 
% the electrodes, headmodel and sourcemodel. 
%
%% visualize the coregistration of electrodes, headmodel, and sourcemodel.
figure(5);

% create colormap to plot parcels in different color
nLabels = length(dkatlas.tissuelabel);
colr = hsv(nLabels); 
vertexcolor = ones(size(dkatlas.pos,1), 3);
for i= 1:length(dkatlas.tissuelabel)
    index = find(dkatlas.tissue==i);
   if ~isempty(index) 
      vertexcolor(index,:) = repmat(colr(i,:),  length(index), 1);
   end   
end

% make the headmodel surface transparent
ft_plot_headmodel(headmodel_eeg, 'edgecolor', 'none','facecolor', 'black'); alpha 0.1
ft_plot_mesh(dkatlas, 'facecolor', 'brain',  'vertexcolor', ...
vertexcolor, 'facealpha', .5);
ft_plot_sens(elec_aligned);
view([0 -90 0])

% Now we can proceed with the computation of the leadfield matrix
cfg = [];
cfg.elec = elec_aligned;            
cfg.channel = dataseg.label;  
cfg.sourcemodel.pos = sourcemodel.pos;              % 2002v source points
cfg.sourcemodel.inside = 1:size(sourcemodel.pos,1); 
cfg.headmodel = headmodel_eeg;                      % volume conduction model
leadfield = ft_prepare_leadfield(cfg);

%% # Source reconstruction and comparison of trials with high and low alpha power
%
% In addition to a forward model, the beamformer needs a sensor-level covariance 
% matrix, or a cross-spectral density matrix. 
% The preliminaries for the cross-spectral density matrix can be obtained with
%
% compute sensor level Fourier spectra, to be used for cross-spectral 
% density computation.
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;
freq           = ft_freqanalysis(cfg, dataseg);

%% do the source reconstruction
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.sourcemodel       = leadfield;
cfg.headmodel         = headmodel_eeg;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori      = 'yes';
cfg.elec              = elec_aligned;
source = ft_sourceanalysis(cfg, freq);
source = ft_sourcedescriptives([], source); % to get the neural-activity-index

%% ## Visualization of the neural-activity-index
% plot the neural activity index (power/noise)

cfg = [];
cfg.parameter    = 'nai';
sourceint = ft_sourceinterpolate(cfg,source,dkatlas);
cfg=[];
sourceint = ft_sourceparcellate(cfg, sourceint, dkatlas);

cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
% cfg.funcolorlim   = [0.0 8];
% cfg.opacitylim    = [3 8];
cfg.opacitymap    = 'rampup';
% cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
figure(4);
ft_sourceplot(cfg, sourceint);
colorbar off
view([-90 30]);
light('Position',[0,-90 30])
material dull
set(gcf,'color','w');

%% ## Creation of a 'pseudo-contrast' based on a median split of the epochs
%
% Typically, in an experimental context, it is useful to visualize activity
% contrasts, e.g., baseline vs. activation intervals, in order to get spatially
% interpretable beamformer results. Although the neural-activity-index intends
% to improve interpretability by normalization with a poor man's approximation
% of the projected noise, and although it takes care of the depth bias of the
% beamformer to some extent, it doesn't usually work well. In order to convince
% ourselves that the beamformer is adequately reconstructing the activity of the
% neural sources, we will resort here to faking an 'experimental' contrast, using
% a median split of the data, where the data are split according to occipital
% alpha power. This requires an estimate of the single epoch alpha power.
% Next, identify the epoch indices for which the alpha power is less/more than
% the median across epochs.
%
% compute sensor level single trial power spectra
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foilim       = [9 11];
cfg.tapsmofrq    = 1;
cfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(cfg, dataseg);
cfg.foilim       = [3 40];
datapowfull           = ft_freqanalysis(cfg, dataseg);

%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind);
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));

%% compute the power spectrum for the median splitted data
cfg              = [];
cfg.trials       = indlow;
datapow_low      = ft_freqdescriptives(cfg, datapowfull);

cfg.trials       = indhigh;
datapow_high     = ft_freqdescriptives(cfg, datapowfull);

%% compute the difference between high and low
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'divide';
powratio      = ft_math(cfg, datapow_high, datapow_low);

%% plot the topography of the difference along with the spectra
cfg        = [];
cfg.layout = lay;
cfg.xlim   = [9.9 10.1];
figure(7);
subplot(1,2,1);ft_topoplotER(cfg, powratio);

cfg         = [];
cfg.channel = {'EEG087', 'EEG088'};
subplot(1,2,2);ft_singleplotER(cfg, datapow_high, datapow_low);

%% ## Source reconstruction of 'low' and 'high' alpha activity epochs
%
% Now we will compute the source reconstructed alpha power again,
% as illustrated above, based on the median split. We will use a common
% filter approach, where we compute the spatial filters based on the
% cross-spectral density averaged across all epochs. 
%
% compute fourier spectra for frequency of interest according to the trial split
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;

cfg.trials = indlow;
freq_low   = ft_freqanalysis(cfg, dataseg);

cfg.trials = indhigh;
freq_high  = ft_freqanalysis(cfg, dataseg);

%% compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.sourcemodel              = leadfield;
cfg.headmodel         = headmodel_eeg;
cfg.elec              = elec_aligned;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, freq);

% use the precomputed filters
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.sourcemodel              = leadfield;
cfg.sourcemodel.filter       = source.avg.filter;
cfg.headmodel         = headmodel_eeg;
cfg.elec              = elec_aligned;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
source_low  = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_low));
source_high = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_high));

cfg           = [];
cfg.operation = 'log10(x1)-log10(x2)';
cfg.parameter = 'pow';
source_ratio  = ft_math(cfg, source_high, source_low);

cfg           = [];
cfg.parameter = 'pow';
sourceint = ft_sourceinterpolate(cfg,source_ratio,dkatlas);
cfg       = [];
sourceint = ft_sourceparcellate(cfg, sourceint, dkatlas);

% We now visualize the log-difference on the cortical sheet.
%
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.colorbar      = 'no';
cfg.funcolormap   = '*RdBu';
figure(8);ft_sourceplot(cfg, sourceint);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);
colorbar off
material dull
set(gcf,'color','w');

% Compare this source reconstruction with the scalp topography generated above. 
% How do the two representations compare?
%
%% # Connectivity analysis and parcellation
%
%% ## Computation of connectivity
% compute connectivity
cfg         = [];
cfg.method  ='coh';
cfg.complex = 'absimag';
source_conn = ft_connectivityanalysis(cfg, source);

% We can now make a, rather uninformative, visualization of the connectome, 
% plotting the full weighted graph, between all pairs of nodes.
%
figure(9);imagesc(source_conn.cohspctrm);

%% ## Parcellation and network analysis
%
% We can now explore the structure in the estimated connectivity matrices using 
% graph theoretic tools. It is not really clear what the effect of the residual 
% spatial leakage of activity is on the estimates of some of these measures, so 
% we would caution for careful interpretations of graph metrics derived from 
% such connectivity matrices, particularly when comparing groups of experimental 
% participants or experimental conditions. Yet, the intention of this tutorial is 
% still to illustrate how such graph theoretic measures can in principle be computed 
% and visualized using fieldtrip. 
cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .1;
network_full = ft_networkanalysis(cfg,source_conn);
%% sourceinterpolate
cfg           = [];
cfg.parameter = 'degrees';
network_int   = ft_sourceinterpolate(cfg,network_full,dkatlas);
cfg           = [];
network_int   = ft_sourceparcellate(cfg, network_int, dkatlas);


% create a fancy mask
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'degrees';
cfg.colorbar      = 'no';
figure(10);ft_sourceplot(cfg, network_int);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);
colorbar off
material dull
set(gcf,'color','w');

% Compare the degree values for the parcellated and the full connectomes. 
% Why are the values different? What determines the maximum value?
%
% Explore and compare both figures a bit more by 3D rotation. 
% Identify the overlap and discrepancies. What could cause this?
%
% Re-compute the node degree based on some other threshold(s), 
% and inspect the effect of threshold on the result.
%
% Re-compute the parcellated connectome using cfg.method = 'max', 
% and inspect the effect of this parameter on the result.
%
%
% Invoke the function and explore the data.
%
%% ## Effect of occipital alpha power on the connectivity results
%
% Compute the connectomes separately on the subsets of trials with low and 
% high occipital alpha power, respectively and inspect the results.
%
%% ## Using other connectivity metrics
%
% Obviously, one can choose from a large amount of different connectivity 
% measures, each of which has its advantages and disadvantages.
%
% Compute the phase locking value between all pairs of dipoles, as well as 
% a parcellated version. Explore the results, and compare them with the 
% imaginary part of coherency.
%
% Compute the envelope correlations using the 'powcorr' method, as well as 
% a parcellated version. Explore the results.
