function test_tutorial_networkanalysis

% MEM 3000mb
% WALLTIME 00:30:00

% TEST test_tutorial_networkanalysis
% TEST ft_networkanalysis

%% read the continuous data and segment into 2 seconds epochs
cfg = [];
cfg.dataset              = dccnpath(fullfile('/home/common/matlab/fieldtrip/data','SubjectRest.ds')); 
cfg.trialdef.triallength = 2;
cfg.trialdef.ntrials     = Inf;

cfg = ft_definetrial(cfg);
cfg.continuous  = 'yes';
cfg.channel     = {'MEG'};
data = ft_preprocessing(cfg);

datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/networkanalysis');
cd(datadir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BYPASS THE INTERACTIVE PART, DECLARE TRIALS 38 91 153 AS BAD
%% make a visual inspection and reject bad trials/sensors
%cfg = [];
%cfg.method  = 'summary';
%cfg.channel = 'MEG';
%cfg.layout  = 'CTF275.lay';
%dataclean = ft_rejectvisual(cfg, data);
%
%% you can check the rejected trial numbers by typing
%trlind = [];
%for i=1:length(dataclean.cfg.artfctdef.summary.artifact)
%  trlind(i) = find(data.sampleinfo(:,1)==dataclean.cfg.artfctdef.summary.artifact(i));
%end;
%disp(trlind);

cfg = [];
cfg.trials = setdiff(1:numel(data.trial), [38 91 153]);
dataclean  = ft_selectdata(cfg, data);

%% downsample the data to speed up component analysis
cfg = [];
cfg.resamplefs = 60;
cfg.detrend    = 'yes';
datads = ft_resampledata(cfg, dataclean);

%% use ICA in order to identify cardiac and blink components
cfg = [];
cfg.method          = 'runica'; 
cfg.runica.maxsteps = 50;
cfg.randomseed      = 0;
comp = ft_componentanalysis(cfg, datads);
%load(fullfile(datadir,'comp.mat'));

%% visualize components

% these were the indices of the bad components that were identified 
% they may be different if you re-run the ICA decomposition
badcomp = [1 3 8 23]; 

cfg = []; 
cfg.channel    = badcomp; 
cfg.layout     = 'CTF275.lay';
cfg.compscale  = 'local';
cfg.continuous = 'yes';
ft_databrowser(cfg, comp);

cfg = [];
cfg.component = badcomp;
dataica = ft_rejectcomponent(cfg, comp);

%% compute the power spectrum
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [1 20];                          
cfg.tapsmofrq    = 2;             
cfg.keeptrials   = 'no';
fft_data = ft_freqanalysis(cfg, dataica);

%% compute the planar transformation
load ctf275_neighb; % loads the neighbourhood structure for the channels

dataicatmp = dataica;
dataicatmp.grad = data.grad;

cfg = [];
cfg.neighbours    = neighbours;
cfg.planarmethod  = 'sincos';
planar = ft_megplanar(cfg, dataicatmp);
clear dataicatmp;

%% compute the power spectrum
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [1 20];                          
cfg.tapsmofrq    = 2;             
cfg.keeptrials   = 'no';
fft_data_planar = ft_freqanalysis(cfg, planar);

%% plot the topography and the spectrum
figure;

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.xlim   = [9 11];
subplot(2,1,1); ft_topoplotER(cfg, fft_data);


cfg = [];
cfg.channel = {'MRO22', 'MRO32', 'MRO33'};
subplot(2,1,2); ft_singleplotER(cfg, fft_data);

%% load the required geometrical information

%% load the template source model, which is in MNI coordinates
%template = load('standard_sourcemodel3d2cm.mat'); 
% load mri % individual mri
% load hdm % individual volume model
% 
% %% compute the source model 
% cfg = [];
% cfg.grid.warpmni   = 'yes';
% cfg.grid.template  = template.sourcemodel;
% cfg.grid.nonlinear = 'yes'; % use non-linear normalization
% cfg.mri            = mri;
% sourcemodel        = ft_prepare_sourcemodel(cfg);

load hdm
load sourcemodel_4k

%% check for the correct alignment of sensors (green) headmodel(transparent) and sourcemodel(blue)
figure;

% make the headmodel surface transparent
ft_plot_vol(hdm, 'edgecolor', 'none');
alpha 0.4           

% add the source model positions and sensors
ft_plot_mesh(ft_convert_units(sourcemodel, 'cm'),'vertexcolor',sourcemodel.sulc);
ft_plot_sens(dataclean.grad);

view([0 -90 0])

%% compute sensor level Fourier spectra
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 2;
cfg.foi        = 10;
freq           = ft_freqanalysis(cfg, dataica);

%% compute the leadfield
cfg             = [];
cfg.grid        = sourcemodel;
cfg.headmodel   = hdm;
cfg.channel     = {'MEG'};
lf = ft_prepare_leadfield(cfg, freq);

%% compute the actual source reconstruction
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.grad              = freq.grad;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.headmodel         = hdm;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori     = 'yes';
source = ft_sourceanalysis(cfg, freq);

%% reduce the source reconstructed data to the dominant orientation
cfg = [];
cfg.projectmom = 'yes';
source_proj = ft_sourcedescriptives(cfg,source);

% and provide the dimension and grid positions of the MRI template positions again
source_proj.dim = template.sourcemodel.dim;
source_proj.pos = template.sourcemodel.pos;

[ftver, ftdir] = ft_version; 
if isunix
   templatefile = [ftdir '/template/anatomy/single_subj_T1.nii'];
elseif ispc
  templatefile =  [ftdir '\template\anatomy\single_subj_T1.nii'];
end

template_mri = ft_read_mri(templatefile);

cfg              = [];
cfg.parameter    = 'nai';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, source_proj, template_mri);

%% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 10];
cfg.opacitylim    = [3 10]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg, source_int);

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'nai';
cfg.location      = [16 -80 -2];
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 10];
cfg.opacitylim    = [3 10]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg,source_int);

%% add the scalp topographies to the figure 

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.xlim   = [9.777650 11.309908];
subplot(2,2,2); ft_topoplotER(cfg,fft_data);
title('axial gradient')

cfg = [];
fft_data_planar_cmb = ft_combineplanar(cfg, fft_data_planar);

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.xlim   = [9.777650 11.309908];
subplot(2,2,3); ft_topoplotER(cfg, fft_data_planar_cmb);
title('planar gradient')

%% compute the power spectrum again but keep the individual trials
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [8 12];                          
cfg.tapsmofrq    = 2;             
cfg.keeptrials   = 'yes';
fft_data = ft_freqanalysis(cfg, dataica);


%% identify the indices of trials with high and low alpha power
tmp     = mean(fft_data.powspctrm,3);           % mean over frequencies between 8-12Hz
ind     = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,ind)<=median(tmp(:,ind)));
indhigh = find(tmp(:,ind)>=median(tmp(:,ind)));

%% compute the planar transformation
cfgneigh          = [];
cfgneigh.feedback = 'no';
cfgneigh.method   = 'template';

cfg = [];
cfg.neighbours    = ft_prepare_neighbours(cfgneigh, dataica);
cfg.planarmethod  = 'sincos';
planar = ft_megplanar(cfg, dataica);

%% compute the power spectrum
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [1 20];                          
cfg.tapsmofrq    = 2;             
cfg.keeptrials   = 'no';

cfg.trials       = indlow; 
tmp = ft_freqanalysis(cfg, planar);
fft_data_planar_low = ft_combineplanar([], tmp);

cfg.trials       = indhigh; 
tmp = ft_freqanalysis(cfg, planar);
fft_data_planar_high = ft_combineplanar([], tmp);

%% and also do the axial representation for comparison purposes
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [1 20];                          
cfg.tapsmofrq    = 2;             
cfg.keeptrials   = 'no';
cfg.trials       = indlow; 
fft_data_low = ft_freqanalysis(cfg,dataica);

cfg.trials       = indhigh; 
fft_data_high = ft_freqanalysis(cfg,dataica);

%% compute the difference between high and low
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract';
diff        = ft_math(cfg, fft_data_high,        fft_data_low);
diff_planar = ft_math(cfg, fft_data_planar_high, fft_data_planar_low);
diff_axial  = ft_math(cfg, fft_data_high,        fft_data_low);

%% plot the topography of the difference along with the spectra
figure;

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.xlim   = [9.777650 11.309908];

subplot(1,3,1); ft_topoplotER(cfg, diff_planar); title ('planar gradient')

subplot(1,3,2); ft_topoplotER(cfg, diff_axial); title ('axial gradient')

cfg = [];
cfg.channel = {'MLO21', 'MRO31'};
subplot(1,3,3); ft_singleplotER(cfg,fft_data_high, fft_data_low);

obj = subplot(1,3,3);
set(obj, 'Position', [.7 .37 .2 .3])
title('')
legend('high alpha','low alpha','Location','northoutside','Orientation','horizontal');

%% compute fourier spectra for frequency of interest
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 2;
cfg.foi        = 10;

cfg.trials     = indlow; 
freq_low = ft_freqanalysis(cfg, dataica);

cfg.trials     = indhigh; 
freq_high = ft_freqanalysis(cfg, dataica);

%% compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.grad              = freq.grad;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.headmodel         = hdm;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '5%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
source = ft_sourceanalysis(cfg, freq);

% use the precomputed filters 
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.grad              = freq.grad;
cfg.method            = 'pcc';
cfg.grid              = lf;
cfg.grid.filter       = source.avg.filter;
cfg.headmodel         = hdm;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '5%';
cfg.pcc.projectnoise  = 'yes';
source_low  = ft_sourceanalysis(cfg, freq_low);
source_high = ft_sourceanalysis(cfg, freq_high);

%% project dipole moments along the dominant orientation
cfg = [];
cfg.projectmom = 'yes';
source_proj_high = ft_sourcedescriptives(cfg,source_high);
source_proj_low  = ft_sourcedescriptives(cfg,source_low);
source_proj_high.dim = template.sourcemodel.dim;
source_proj_high.pos = template.sourcemodel.pos;
source_proj_low.dim  = template.sourcemodel.dim;
source_proj_low.pos  = template.sourcemodel.pos;

%% interpolate the results onto an anatomical brain template
[ftver, ftdir] = ft_version; 
if isunix
   templatefile = [ftdir '/template/anatomy/single_subj_T1.nii'];
elseif ispc
  templatefile =  [ftdir '\template\anatomy\single_subj_T1.nii'];
end

template_mri = ft_read_mri(templatefile);

cfg              = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int_high  = ft_sourceinterpolate(cfg, source_proj_high, template_mri);
source_int_low   = ft_sourceinterpolate(cfg, source_proj_low,  template_mri);

cfg  = [];
cfg.operation = '(x1-x2)/x2';
cfg.parameter = 'pow';
source_int = ft_math(cfg, source_int_high, source_int_low);

cfg  = [];
cfg.operation = '(x1-x2)/x2';
cfg.parameter = 'pow';
source_int = ft_math(cfg, source_int_high, source_int_low);

% up to 50 percent of maximum
source_int.mask = source_int.pow > max(source_int.pow(:))*.5;

% copy the anatomy
source_int.anatomy = source_int_high.anatomy;

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-1 1];
cfg.location      = [16 -78 38];
cfg.funcolormap   = 'jet';
ft_sourceplot(cfg, source_int);

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.xlim   = [9.777650 11.309908];
subplot(2,2,2); ft_topoplotER(cfg, diff_planar);

%% reduce memory demands and compute connectivity

% compute the sparse representation
source_sparse = ft_source2sparse(source_proj);

% then compute connectivity
cfg=[];
cfg.method  ='coh';
cfg.complex = 'absimag';
source_conn = ft_connectivityanalysis(cfg, source_sparse);

%% reassign the connectivity
cohspctrm_full = nan(size(source_proj.pos,1));
for i=1:size(source_conn.pos)
  pos1 = source_conn.pos(i,1:3);
  pos2 = source_conn.pos(i,4:6);
  ind1 = find(source_proj.pos(:,1)==pos1(1) & source_proj.pos(:,2)==pos1(2) & source_proj.pos(:,3)==pos1(3));
  ind2 = find(source_proj.pos(:,1)==pos2(1) & source_proj.pos(:,2)==pos2(2) & source_proj.pos(:,3)==pos2(3));
  cohspctrm_full(ind1,ind2) = source_conn.cohspctrm(i);
end

% then go to the 'full' representation again
source_conn_full = [];
source_conn_full.pos       = source_proj.pos;
source_conn_full.dim       = source_proj.dim;
source_conn_full.inside    = source_proj.inside;
source_conn_full.cohspctrm = cohspctrm_full;
source_conn_full.dimord    = 'pos_pos';

%% compute graph metric
cfg = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .1;
network = ft_networkanalysis(cfg,source_conn_full);

%% interpolate on anatomical mri and plot the result
network.pos     = template.sourcemodel.pos;
network.dim     = template.sourcemodel.dim;
network.inside  = template.sourcemodel.inside;

cfg              = [];
cfg.parameter    = 'degrees';
cfg.interpmethod = 'nearest';
network_int  = ft_sourceinterpolate(cfg, network, template_mri);

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

node = zeros(372,6);
node(:,1:3) = sourcemodel.pos(sourcemodel.inside,:);
node(:,4)   = 4; 
node(:,5)   = network.degrees(network.inside);
node(:,6)   = 0;
node(:,1:3) = node(:,1:3)*10;
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


