function failed_bug1049

% MEM 2000mb
% WALLTIME 00:10:00

% TEST ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis

% this function creates a set of source-structures to be used for testing

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];
ft_default.feedback = 'no';

% get volume conductor model
volname = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151/Subject01.ds/default.hdm');
vol     = ft_read_vol(volname);

% get data + sensor info
dataname = dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat');
load(dataname);

% create 3D grid
cfg      = [];
cfg.grad = data.grad;
cfg.vol  = vol;
cfg.channel = 'MEG';
cfg.grid.resolution = 1.5;
grid = ft_prepare_leadfield(cfg);

% create 2D grid
[pnt, tri] = icosahedron162;
pnt   = pnt*(vol.orig.MEG_Sphere.RADIUS-1.5);
shift = [vol.orig.MEG_Sphere.ORIGIN_X vol.orig.MEG_Sphere.ORIGIN_Y vol.orig.MEG_Sphere.ORIGIN_Z];
pnt   = pnt+repmat(shift,[size(pnt,1) 1]);
grid2.pnt = pnt;
grid2.tri = tri;
grid2.inside = 1:size(grid2.pnt,1);
grid2.outside = [];

cfg      = [];
cfg.grad = data.grad;
cfg.vol  = vol;
cfg.grid = grid2;
cfg.channel = 'MEG';
grid2 = ft_prepare_leadfield(cfg);

% create timelock structure with covariance for lcmv beamforming and
% minimumnormestimate
cfg  = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
cfg.channel    = 'MEG';
tlck = ft_timelockanalysis(cfg, data);

% create freq structure for dics beamforming and pcc beamforming
cfg  = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.tapsmofrq = 4;
cfg.foilim = [0 20];
cfg.channel = 'MEG';
freq = ft_freqanalysis(cfg, data);

outputdir = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/source/meg'));

% do LCMV beamforming
cfg            = [];
cfg.method     = 'lcmv';
cfg.lcmv.keepleadfield = 'yes';
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.keepcov       = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.grid       = grid;
cfg.vol        = vol;
cfg.outputfile = fullfile(outputdir, 'ctf151_lcmv3d_avg');
sourcelcmv3d1  = ft_sourceanalysis(cfg, tlck);
cfg.grid       = grid2;
cfg.outputfile = 'ctf151_lcmv2d_avg';
sourcelcmv2d1  = ft_sourceanalysis(cfg, tlck);

cfg.rawtrial    = 'yes';
cfg.grid        = grid;
cfg.outputfile  = fullfile(outputdir, 'ctf151_lcmv3d_trial');
cfg.grid.filter = sourcelcmv3d1.avg.filter;
ft_sourceanalysis(cfg, tlck);
cfg.grid        = grid2;
cfg.outputfile  = fullfile(outputdir, 'ctf151_lcmv2d_trial');
cfg.grid.filter = sourcelcmv2d1.avg.filter;
ft_sourceanalysis(cfg, tlck);

% do MNE
cfg = [];
cfg.method   = 'mne';
cfg.mne.keepleadfield = 'yes';
cfg.mne.keepfilter = 'yes';
cfg.mne.lambda     = 1e4;
cfg.vol  = vol;
cfg.grid = grid;
cfg.outputfile = fullfile(outputdir, 'ctf151_mne3d');
sourcemne3d1 = ft_sourceanalysis(cfg, tlck);
cfg.grid = grid2;
cfg.outputfile = fullfile(outputdir, 'ctf151_mne2d');
sourcemne2d1 = ft_sourceanalysis(cfg, tlck);

cfg.rawtrial    = 'yes';
cfg.grid        = grid;
cfg.grid.filter = sourcemne3d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_mne3d_trial');
ft_sourceanalysis(cfg, tlck);
cfg.grid        = grid2;
cfg.grid.filter = sourcemne2d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_mne2d_trial');
ft_sourceanalysis(cfg, tlck);

% do DICS
cfg = [];
cfg.method = 'dics';
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.lambda        = '5%';
cfg.frequency = 10;
cfg.vol  = vol;
cfg.grid = grid;
cfg.outputfile = fullfile(outputdir, 'ctf151_dics3d_avg');
sourcedics3d1 = ft_sourceanalysis(cfg, freq);
cfg.grid = grid2;
cfg.outputfile = fullfile(outputdir, 'ctf151_dics2d_avg');
sourcedics2d1 = ft_sourceanalysis(cfg, freq);

cfg.rawtrial    = 'yes';
cfg.grid        = grid;
cfg.grid.filter = sourcedics3d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_dics3d_trial');
ft_sourceanalysis(cfg, freq);
cfg.grid        = grid2;
cfg.grid.filter = sourcedics2d1.avg.filter;
cfg.outputfile  = fullfile(outputdir, 'ctf151_dics2d_trial');
ft_sourceanalysis(cfg, freq);


% do PCC
cfg = [];
cfg.method = 'pcc';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.keepleadfield = 'yes';
cfg.pcc.keepcsd       = 'yes';
cfg.pcc.keepmom       = 'yes';
cfg.pcc.lambda        = '5%';
cfg.frequency = 10;
cfg.vol  = vol;
cfg.grid = grid;
cfg.outputfile = fullfile(outputdir, 'ctf151_pcc3d');
ft_sourceanalysis(cfg, freq);
cfg.grid = grid2;
cfg.outputfile = fullfile(outputdir, 'ctf151_pcc2d');
ft_sourceanalysis(cfg, freq);


