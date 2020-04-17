function test_pull1377
% MEM 2gb
% WALLTIME 00:10:00

% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis 

% this function creates a set of source-structures to be used for testing

% get volume conductor model
% for MEG, localsphere
volname = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151/Subject01.ds/default.hdm');
% figure out what model this is, preferably local spheres
vol     = ft_read_headmodel(volname);

% for EEG, singlesphere

% get MEG data + sensor info
dataname = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');
load(dataname);

% get EEG data + channel info

% create 3D grid
% for MEG
cfg      = [];
cfg.grad = data.grad;
cfg.headmodel = vol;
cfg.channel = 'MEG';
cfg.sourcemodel.resolution = 1.5;
grid = ft_prepare_leadfield(cfg);

% create 2D grid
[pnt, tri] = mesh_sphere(162);
pnt   = pnt*(vol.orig.MEG_Sphere.RADIUS-1.5);
shift = [vol.orig.MEG_Sphere.ORIGIN_X vol.orig.MEG_Sphere.ORIGIN_Y vol.orig.MEG_Sphere.ORIGIN_Z];
pnt   = pnt+repmat(shift,[size(pnt,1) 1]);
grid2.pnt = pnt;
grid2.tri = tri;
grid2.inside = 1:size(grid2.pnt,1);
grid2.outside = [];

cfg      = [];
cfg.grad = data.grad;
cfg.headmodel = vol;
cfg.sourcemodel = grid2;
cfg.channel = 'MEG';
grid2 = ft_prepare_leadfield(cfg);

% for EEG

% load externally created leadfields
% for MEG

% for EEG

% create timelock structure with covariance for lcmv beamforming and
% minimumnormestimate
% for MEG
cfg  = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
cfg.channel    = 'MEG';
MEG_tlck = ft_timelockanalysis(cfg, data);

% for EEG 

% create freq structure for dics beamforming and pcc beamforming
% for MEG,
cfg  = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.tapsmofrq = 4;
cfg.foilim = [0 20];
cfg.channel = 'MEG';
MEG_freq = ft_freqanalysis(cfg, data);

% for EEG

outputdir = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/source/meg'));

% inverse solutions

% do LCMV beamforming
cfg            = [];
cfg.method     = 'lcmv';
cfg.lcmv.keepleadfield = 'yes';
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.keepcov       = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.sourcemodel       = grid;
cfg.headmodel  = vol;
cfg.outputfile = fullfile(outputdir, 'ctf151_lcmv3d_avg');
sourcelcmv3d1  = ft_sourceanalysis(cfg, tlck);
cfg.sourcemodel       = grid2;
cfg.outputfile = 'ctf151_lcmv2d_avg';
sourcelcmv2d1  = ft_sourceanalysis(cfg, MEG_tlck);

cfg.rawtrial    = 'yes';
cfg.sourcemodel        = grid;
cfg.outputfile  = fullfile(outputdir, 'ctf151_lcmv3d_trial');
cfg.sourcemodel.filter = sourcelcmv3d1.avg.filter;
ft_sourceanalysis(cfg, tlck);
cfg.sourcemodel        = grid2;
cfg.outputfile  = fullfile(outputdir, 'ctf151_lcmv2d_trial');
cfg.sourcemodel.filter = sourcelcmv2d1.avg.filter;
ft_sourceanalysis(cfg, MEG_tlck);

% do MNE
cfg = [];
cfg.method   = 'mne';
cfg.mne.keepleadfield = 'yes';
cfg.mne.keepfilter = 'yes';
cfg.mne.lambda     = 1e4;
cfg.headmodel = vol;
cfg.sourcemodel = grid;
cfg.outputfile = fullfile(outputdir, 'ctf151_mne3d');
sourcemne3d1 = ft_sourceanalysis(cfg, tlck);
cfg.sourcemodel = grid2;
cfg.outputfile = fullfile(outputdir, 'ctf151_mne2d');
sourcemne2d1 = ft_sourceanalysis(cfg, MEG_tlck);

cfg.rawtrial    = 'yes';
cfg.sourcemodel        = grid;
cfg.sourcemodel.filter = sourcemne3d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_mne3d_trial');
ft_sourceanalysis(cfg, tlck);
cfg.sourcemodel        = grid2;
cfg.sourcemodel.filter = sourcemne2d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_mne2d_trial');
ft_sourceanalysis(cfg, MEG_tlck);

% do DICS
cfg = [];
cfg.method = 'dics';
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.lambda        = '5%';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.sourcemodel = grid;
cfg.outputfile = fullfile(outputdir, 'ctf151_dics3d_avg');
sourcedics3d1 = ft_sourceanalysis(cfg, MEG_freq);
cfg.sourcemodel = grid2;
cfg.outputfile = fullfile(outputdir, 'ctf151_dics2d_avg');
sourcedics2d1 = ft_sourceanalysis(cfg, MEG_freq);

cfg.rawtrial    = 'yes';
cfg.sourcemodel        = grid;
cfg.sourcemodel.filter = sourcedics3d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_dics3d_trial');
ft_sourceanalysis(cfg, MEG_freq);
cfg.sourcemodel        = grid2;
cfg.sourcemodel.filter = sourcedics2d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_dics2d_trial');
ft_sourceanalysis(cfg, MEG_freq);


% do PCC
cfg = [];
cfg.method = 'pcc';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.keepleadfield = 'yes';
cfg.pcc.keepcsd       = 'yes';
cfg.pcc.keepmom       = 'yes';
cfg.pcc.lambda        = '5%';
cfg.frequency = 10;
cfg.headmodel = vol;
cfg.sourcemodel = grid;
cfg.outputfile = fullfile(outputdir, 'ctf151_pcc3d');
ft_sourceanalysis(cfg, freq);
cfg.sourcemodel = grid2;
cfg.outputfile = fullfile(outputdir, 'ctf151_pcc2d');
ft_sourceanalysis(cfg, MEG_freq);

% do dipolefit
% for MEG
cfg = [];
cfg.numdipoles    = 1;                              %number of expected sources
cfg.headmodel     = headmodelr;                     %the head model
cfg.grid          = leadfield;                      %the (precomputed) leadfield
cfg.nonlinear     = 'no';                           %only dipole scan
cfg.grad          = grad;                           %the sensor model
cfg.latency       = 0.025;                          %the latency of interest
%dipfit_meg    = ft_dipolefitting(cfg,MEG_avg);
ft_dipolefitting(cfg,MEG_tlck); % we do not want to save the output I guess

% for EEG
cfg = [];
cfg.numdipoles    =  1;                             
cfg.headmodel     = headmodel;                  
cfg.grid          = leadfield;    
cfg.nonlinear     = 'no';                          
cfg.elec          = elec;                           
cfg.latency       = 0.025;                         
ft_dipolefitting(cfg,EEG_tlck);

