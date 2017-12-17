function test_bug2511

% WALLTIME 00:20:00
% MEM 8gb

% TEST ft_sourceplot ft_read_headshape

t1 = ft_read_mri(dccnpath('/home/common/matlab/spm8/canonical/single_subj_T1.nii'));
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2511.mat'));

% loading the data takes quite some time, as it is 4.7GB on disk, which is
% a bit silly, because only 1 variable is used, and the whole
% XXX.cfg.previous...... is not needed. JM reduced the file size on
% 20140629

cfg                  = [];
cfg.parameter        = 'avg.pow';
cfg.interpmethod    = 'sphere_avg';
sourceDiffIntAVG     = ft_sourceinterpolate(cfg, AverageSource ,t1);

cfgS                     = [];
cfgS.funparameter        = 'avg.pow';
cfgS.funcolormap         = 'jet';
cfgS.opacitylim          = [0 1];
cfgS.opacitymap          = 'rampup';
cfgS.method              = 'surface';
cfgS.projmethod          = 'nearest';
cfgS.interactive         = 'yes';
% cfgS.surffile            = 'surface_l4_both.mat'; % this was causing the problem, the file was removed from FT
cfgS.surffile            = 'surface_white_both.mat'; % this is the new surface that should be used instead
cfgS.surfdownsample      = 10;
cfgS.funcolorlim         = 'maxabs';
cfgS.interactive         = 'yes';
ft_sourceplot(cfgS, sourceDiffIntAVG);

%% error occurs with
% addpath(dccnpath('/home/common/matlab/fieldtrip/'));
%% but not with
% addpath(dccnpath('/home/common/matlab/fieldtrip-20131231/'));

%% Error description
%
% the input is volume data with dimensions [91 109 91]
% scaling anatomy to [0 1]
% no masking parameter
% The source data is defined on a 3D grid, interpolation to a surface mesh will be performed
% Warning: could not determine filetype of surface_l4_both.mat
% > In fileio/private/ft_warning at 158
%   In ft_filetype at 1189
%   In ft_read_headshape at 110
%   In ft_sourceplot at 890
% Error using ft_read_headshape (line 888)
% unknown fileformat "unknown" for head shape information
%
% Error in ft_sourceplot (line 890)
%     surf  = ft_read_headshape(cfg.surffile);
