function test_tutorial_networkanalysis_new

% MEM 3000mb
% WALLTIME 00:30:00

% TEST test_tutorial_networkanalysis
% TEST ft_networkanalysis

%% read the continuous data and segment into 2 seconds epochs, with 50% overlap
cfg            = [];
cfg.dataset    = dccnpath(fullfile('/home/common/matlab/fieldtrip/data','SubjectRest.ds')); 
cfg.continuous = 'yes';
cfg.channel    = {'MEG'};
data = ft_preprocessing(cfg);

cfg         = [];
cfg.length  = 2;
cfg.overlap = 0.5;
data        = ft_redefinetrial(cfg, data);

cfg        = [];
cfg.demean = 'yes';
cfg.trials = 1:(numel(data.trial)-6);
data       = ft_preprocessing(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BYPASS THE INTERACTIVE PART, DECLARE TRIALS 38 91 153 AS BAD
%% make a visual inspection and reject bad trials/sensors
% cfg = [];
% cfg.method  = 'summary';
% cfg.channel = 'MEG';
% cfg.layout  = 'CTF275.lay';
% dataclean = ft_rejectvisual(cfg, data);
%
%% you can check the rejected trial numbers by typing
% trlind = [];
% for i=1:length(dataclean.cfg.artfctdef.summary.artifact)
%  trlind(i) = find(data.sampleinfo(:,1)==dataclean.cfg.artfctdef.summary.artifact(i));
% end;
% disp(trlind);

badtrials  = [18 19 21 72 73 74 75 76 93 94 109 110 126 127 128 140 172 173 179 180 181 182 196 197 198 227 228 233 243 244 250 251 265 266 286];
cfg        = [];
cfg.trials = setdiff(1:numel(data.trial), badtrials);
dataclean  = ft_selectdata(cfg, data);

%% downsample the data to speed up component analysis
dataclean.time(1:end) = dataclean.time(1);

cfg            = [];
cfg.resamplefs = 100;
cfg.detrend    = 'yes';
datads         = ft_resampledata(cfg, dataclean);

%% use ICA in order to identify cardiac and blink components
cfg                 = [];
cfg.method          = 'runica'; 
cfg.runica.maxsteps = 50;
cfg.randomseed      = 0;
comp                = ft_componentanalysis(cfg, datads);

%% visualize components

% these were the indices of the bad comp**[[reference:ft_definetrial|ft_definetrial]]** and onents that were identified 
% they may be different if you re-run the ICA decomposition
badcomp = [2 3 7 16]; 
 
cfg = []; 
cfg.channel    = badcomp; 
cfg.layout     = 'CTF275_helmet.mat';
cfg.compscale  = 'local';
cfg.continuous = 'yes';
ft_databrowser(cfg, comp);

cfg           = [];
cfg.component = badcomp;
dataica       = ft_rejectcomponent(cfg, comp);

%% compute the power spectrum
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'no';
datapow          = ft_freqanalysis(cfg, dataica);

%% compute the planar transformation, this is not really necessary, but instructive anyhow
load ctf275_neighb; % loads the neighbourhood structure for the channels

dataicatmp      = dataica;
dataicatmp.grad = datads.grad;

cfg               = [];
cfg.neighbours    = neighbours;
cfg.planarmethod  = 'sincos';
planar            = ft_megplanar(cfg, dataicatmp);
clear dataicatmp;

%% compute the power spectrum
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'no';
datapow_planar   = ft_freqanalysis(cfg, planar);

%% plot the topography and the spectrum
figure;

cfg        = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.xlim   = [9 11];
subplot(2,2,1); ft_topoplotER(cfg, datapow);
subplot(2,2,2); ft_topoplotER(cfg, datapow_planar);

cfg         = [];
cfg.channel = {'MRO22', 'MRO32', 'MRO33'};
subplot(2,2,3); ft_singleplotER(cfg, datapow);


%% load the required geometrical information
load hdm
load sourcemodel_4k

%% visualize the coregistration of sensors, headmodel, and sourcemodel.
figure;

% make the headmodel surface transparent
ft_plot_vol(hdm, 'edgecolor', 'none'); alpha 0.4           
ft_plot_mesh(ft_convert_units(sourcemodel, 'cm'),'vertexcolor',sourcemodel.sulc);
ft_plot_sens(dataclean.grad);
view([0 -90 0])

%% compute the leadfield
cfg             = [];
cfg.grid        = sourcemodel;
cfg.headmodel   = hdm;
cfg.channel     = {'MEG'};
lf              = ft_prepare_leadfield(cfg, dataica);

%% compute sensor level Fourier spectra
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;
freq           = ft_freqanalysis(cfg, dataica);


%% compute the actual source reconstruction
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.headmodel         = hdm;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, freq);
source = ft_sourcedescriptives([], source);

%% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 8];
cfg.opacitylim    = [3 8]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source);
view([-90 30]);
light;

%% compute the power spectrum again but keep the individual trials
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [0 30];                          
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'yes';
datapow = ft_freqanalysis(cfg, dataica);


%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind);    
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));

%% compute the power spectrum for the median splitted data
cfg              = [];
cfg.trials       = indlow; 
datapow_low      = ft_freqdescriptives(cfg, datapow);

cfg.trials       = indhigh; 
datapow_high     = ft_freqdescriptives(cfg, datapow);

%% compute the difference between high and low
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'divide';
powratio      = ft_math(cfg, datapow_high, datapow_low);

%% plot the topography of the difference along with the spectra
cfg        = [];
cfg.layout = 'CTF275_helmet.mat';
cfg.xlim   = [9.9 10.1];
figure; ft_topoplotER(cfg, powratio);

cfg         = [];
cfg.channel = {'MRO33'};
figure; ft_singleplotER(cfg, datapow_high, datapow_low);


%% compute fourier spectra for frequency of interest according to the trial split
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;

cfg.trials = indlow; 
freq_low   = ft_freqanalysis(cfg, dataica);

cfg.trials = indhigh; 
freq_high  = ft_freqanalysis(cfg, dataica);

%% compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.headmodel         = hdm;
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
cfg.grid              = lf;
cfg.grid.filter       = source.avg.filter;
cfg.headmodel         = hdm;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
source_low  = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_low));
source_high = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_high));

cfg           = [];
cfg.operation = 'log10(x1)-log10(x2)';
cfg.parameter = 'pow';
source_ratio  = ft_math(cfg, source_high, source_low);

% create a fancy mask
source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2; 

cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_ratio);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);

% compute connectivity
cfg         = [];
cfg.method  ='coh';
cfg.complex = 'absimag';
source_conn = ft_connectivityanalysis(cfg, source);

% parcellate
load atlas_MSMAll_4k
atlas.pos = source_conn.pos; % otherwise the parcellation won't work

cfg = [];
cfg.parcellation = 'parcellation';
cfg.parameter    = 'cohspctrm';
parc_conn = ft_sourceparcellate(cfg, source_conn, atlas);


%% compute some graph metric
cfg = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .1;
network_full = ft_networkanalysis(cfg,source_conn);
network_parc = ft_networkanalysis(cfg,parc_conn);

%%
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'avg.degrees';
cfg.funcolormap   = 'jet';
cfg.location      = 'max';
cfg.funparameter  = 'degrees';
ft_sourceplot(cfg, network_int);

% we increase the threshold here to highlight the dominant links
edge = source_conn.cohspctrm >.15;
dlmwrite('edge.edge',edge,'\t');

load standard_sourcemodel3d2cm.mat

node = zeros(372,6)
node(:,1:3) = sourcemodel.pos(sourcemodel.inside,:)
node(:,4)   = 4; 
node(:,5)   = network.degrees(network.inside);
node(:,6)   = 0;
node(:,1:3) = node(:,1:3)*10
dlmwrite('node.node',node,' ');

BrainNet_MapCfg('/yourpath/mesh.nv','/yourpath/node.node','/yourpath/edge.edge');
view([0 -90 0])

cfg = [];
cfg.method = 'plv';
source_conn = ft_connectivityanalysis(cfg, source_sparse);

%% then compute graph metric
source_conn.dim     = source_sparse.dim;
source_conn.outside = source_proj.outside;

% then go to the 'full' representation again
source_conn_full = ft_source2full(source_conn);
source_conn_full.dimord='pos_pos';

%%
fn=fieldnames(source_conn_full);
cfg = [];
cfg.method = 'degrees';
cfg.parameter = fn{4};
cfg.threshold = .5;
deg = ft_networkanalysis(cfg,source_conn_full);


