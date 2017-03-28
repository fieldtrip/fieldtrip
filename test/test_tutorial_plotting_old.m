function test_tutorial_plotting

% WALLTIME 00:15:00
% MEM 4gb

% TEST ft_multiplotER ft_singleplotER ft_topoplotER ft_singleplotTFR ft_multiplotTRF ft_megplanar ft_combineplanar ft_volumereslice

% this tutorial reflects the plotting tutorial up to March 2017
% around that time we did a data visualization workshop at the donders for which I created an updated tutorial

% see http://www.fieldtriptoolbox.org/tutorial/plotting
% this testscript corresponds to the version on the wiki at 23 December 2012

% use the tutorial dataset from home/common

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/plotting/avgFC.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/plotting/GA_FC.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/plotting/TFRhann.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/plotting/statERF.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/plotting/statTFR.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/plotting/sourceDiff.mat'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Singleplot Functions %%%%%%%%%
% ft_singleplotER
cfg = [];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1e-13 3e-13];
cfg.channel = 'MLC24';
clf;
ft_singleplotER(cfg,avgFC);

% MATLAB toolbox alternative of ft_singleplotER
selected_data = avgFC.avg(9,241:601); %MLC24 is the 9th channel, -0.2 to 1.0 is sample 241 to 601
time = avgFC.time(241:601);
figure;
plot(time, selected_data)
xlim([-0.2 1.0])
ylim([-1e-13 3e-13])


% ft_singleplotTFR
cfg = [];
cfg.baseline = [-0.5 -0.1]; 
cfg.baselinetype = 'absolute'; 	
cfg.zlim = [-1.5e-27 1.5e-27];	
cfg.channel       = 'MLC*';
cfg.renderer      = 'painters'; % might be needed if colorbar does not print properly
figure;ft_singleplotTFR(cfg, TFRhann); colorbar;

cfg.channel       = 'MRC*';
figure;ft_singleplotTFR(cfg, TFRhann); colorbar;


%%%%%%%%% Topoplot Functions %%%%%%%%%
cfg = [];                            
cfg.xlim = [0.3 0.5];                
cfg.zlim = [0 6e-14];                
cfg.layout = 'CTF151.lay';            
figure; ft_topoplotER(cfg,GA_FC); colorbar;    

cfg = [];
cfg.xlim = [0.9 1.3];                
cfg.ylim = [15 20];                  
cfg.zlim = [-1e-27 1e-27];           
cfg.baseline = [-0.5 -0.1];          
cfg.baselinetype = 'absolute';
cfg.layout = 'CTF151.lay';
figure; ft_topoplotTFR(cfg,TFRhann)

% multiple topoplots in one figure
cfg.xlim = -0.4:0.2:1.4;
cfg.comment = 'xlim';
cfg.commentpos = 'title';
figure; ft_topoplotTFR(cfg,TFRhann)

% data selection options prior to calling topoplot
% not specific to topoplot
cfg = [];
cfg.xlim = [0.9 1.3];
cfg.ylim = [15 20];
cfg.zlim = [-1e-27 1e-27];
cfg.baseline = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.layout = 'CTF151.lay';

% specific to topoplot - option 1 (regular colours)
cfg.gridscale = 300;                  
cfg.style = 'straight';               
cfg.marker = 'labels';                
figure; ft_topoplotTFR(cfg,TFRhann)

% specific to topoplot - option 2 (gray scale)
cfg.gridscale = 300;                
cfg.contournum = 10;                
cfg.colormap = gray(10);            
figure; ft_topoplotTFR(cfg,TFRhann)

% specific to topoplot - option 3 (spring colormap)
cfg.gridscale = 300;
cfg.contournum = 4;
cfg.colormap = spring(4);
cfg.markersymbol = '.';
cfg.markersize = 12;
cfg.markercolor = [0 0.69 0.94];
figure; ft_topoplotTFR(cfg,TFRhann)

%%%%%%%%% Plotting with interactive mode / multiplots  %%%%%%%%%

cfg = [];
cfg.baseline = [-0.5 -0.1];
cfg.zlim = [-3e-27 3e-27];
cfg.baselinetype = 'absolute';
cfg.layout = 'CTF151.lay';
cfg.interactive = 'yes';
figure; ft_multiplotTFR(cfg,TFRhann)

cfg = [];
cfg.xlim = [-0.2 1.0];
cfg.ylim = [-1e-13 3e-13];
cfg.layout = 'CTF151.lay';
cfg.interactive = 'yes';
figure; ft_multiplotER(cfg,avgFC);


%%%%%%%%% Connectivity plots  %%%%%%%%%

cfg = [];
cfg.covariance = 'yes';
tmp1 = ft_timelockanalysis(cfg, avgFC);

cfg = [];
cfg.method = 'corr';
tmp2 = ft_connectivityanalysis(cfg, tmp1);

cohFC      = tmp2;
cohFC.coh  = tmp2.corr;
cohFC.freq = 1;
cohFC = rmfield(cohFC, 'coh');

cfg = [];
cfg.foi = 1;
cfg.layout = 'CTF151.lay';
figure; ft_topoplotCC(cfg, cohFC);

%%%%%%%%% Cluster plots  %%%%%%%%%

% timelocked data 
cfg = [];
cfg.zlim = [-6 6]; %Tvalues
cfg.alpha = 0.05;
cfg.layout = 'CTF151.lay';
figure; ft_clusterplot(cfg,statERF);

% freq data 
cfg = [];
cfg.zlim = [-5 5];
cfg.alpha = 0.05;
cfg.layout = 'CTF151.lay';
ft_clusterplot(cfg,statTFR);

%%%%%%%% beamforming plotting  %%%%%%%

% (1) multiple axial slices

%% load MRI and interpolate functional source data to MRI
mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/plotting/Subject01.mri'));  
mri = ft_volumereslice([], mri);

cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'pow';
sourceDiffInt  = ft_sourceinterpolate(cfg, sourceDiff , mri);

%% plot multiple 2D axial slices
cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2]; 
cfg.opacitymap    = 'rampup';  
ft_sourceplot(cfg, sourceDiffInt);

%% (2) ortho plot
% volume normalise
cfg = [];
cfg.coordsys      = 'ctf';
cfg.nonlinear     = 'no';
sourceDiffIntNorm = ft_volumenormalise(cfg, sourceDiffInt);

% plot ortho
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2]; 
cfg.opacitymap    = 'rampup';  
figure; ft_sourceplot(cfg, sourceDiffIntNorm);

%% (3) MNI surface
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolorlim    = [0.0 1.2];
cfg.funcolormap    = 'jet';
cfg.opacitylim     = [0.0 1.2]; 
cfg.opacitymap     = 'rampup';  
cfg.projmethod     = 'nearest'; 
cfg.surffile      = 'surface_white_both.mat';
cfg.surfdownsample = 10; 
ft_sourceplot(cfg, sourceDiffIntNorm);
view ([90 0])
