function test_tutorial_networkanalysis(datadirs)

% WALLTIME 00:30:00
% MEM 7gb

% TEST ft_networkanalysis

if nargin==0,
    datadirs{1} = dccnpath('/home/common/matlab/fieldtrip/data');
    datadirs{2} = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/networkanalysis');
end

%% read the continuous data and segment into 2 seconds epochs, with 50% overlap
cfg            = [];
cfg.dataset    = fullfile(datadirs{1},'SubjectRest.ds'); 
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

datadir = datadirs{2};
cd(datadir);

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
%cfg                 = [];
%cfg.method          = 'runica'; 
%cfg.runica.maxsteps = 50;
%cfg.randomseed      = 0;
%comp                = ft_componentanalysis(cfg, datads);
load(fullfile(datadir,'comp.mat'));

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
subplot(2,2,2); ft_topoplotER(cfg, ft_combineplanar([], datapow_planar));

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

figure;imagesc(source_conn.cohspctrm);

% parcellate
load atlas_MMP1.0_4k.mat
atlas.pos = source_conn.pos; % otherwise the parcellation won't work

cfg = [];
cfg.parcellation = 'parcellation';
cfg.parameter    = 'cohspctrm';
parc_conn = ft_sourceparcellate(cfg, source_conn, atlas);

figure;imagesc(parc_conn.cohspctrm);


%% compute node degree
cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .1;
network_full = ft_networkanalysis(cfg,source_conn);
network_parc = ft_networkanalysis(cfg,parc_conn);

%% visualize
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'degrees';
cfg.funcolormap   = 'jet';
cfg.location      = 'max';
cfg.funparameter  = 'degrees';
figure; ft_sourceplot(cfg, network_full);
view([-150 30]);

figure; ft_sourceplot(cfg, network_parc);
view([-150 30]);

%% now inspect some (averages of) columns of the connectivity matrix
load atlas_MMP1.0_4k.mat;
load sourcemodel_4k_inflated;


%% alpha power dependent connectivity
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind);    
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
[sortpow, ix] = sort(tmp(:,chanind));

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

tmpcfg      = [];
tmpcfg.trials  = ix(1:60);
source_low  = ft_sourcedescriptives([], ft_sourceanalysis(cfg, ft_selectdata(tmpcfg, freq)));
tmpcfg.trials  = ix(end-60:end);
source_high = ft_sourcedescriptives([], ft_sourceanalysis(cfg, ft_selectdata(tmpcfg, freq)));


cfg              = [];
cfg.method       = 'coh';
%cfg.complex      = 'absimag';
source_conn_low  = ft_connectivityanalysis(cfg, source_low);
source_conn_high = ft_connectivityanalysis(cfg, source_high);

atlas.pos = source_conn_low.pos;

cfg = [];
cfg.parcellation = 'parcellation';
cfg.parameter    = 'cohspctrm';
parc_conn_low    = ft_sourceparcellate(cfg, source_conn_low,  atlas);
parc_conn_high   = ft_sourceparcellate(cfg, source_conn_high, atlas); 

cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .1;
network_full_high = ft_networkanalysis(cfg,source_conn_high);
network_full_low  = ft_networkanalysis(cfg,source_conn_low);
network_parc_high = ft_networkanalysis(cfg,parc_conn_high);
network_parc_low  = ft_networkanalysis(cfg,parc_conn_low);

cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'degrees';
cfg.funcolormap   = 'jet';
cfg.location      = 'max';
cfg.funparameter  = 'degrees';
figure; ft_sourceplot(cfg, network_full_high);
view([-150 30]);
figure; ft_sourceplot(cfg, network_full_low);
view([-150 30]);

figure; ft_sourceplot(cfg, network_parc_high);
view([-150 30]);
figure; ft_sourceplot(cfg, network_parc_low);
view([-150 30]);
