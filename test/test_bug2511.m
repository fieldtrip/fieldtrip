function test_bug2511

% WALLTIME 00:20:00
% MEM 8gb

% TEST test_bug2511
% TEST ft_sourceplot ft_read_headshape

t1                   = ft_read_mri('/home/common/matlab/spm8/canonical/single_subj_T1.nii');
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
cfgS.surffile            = 'surface_l4_both.mat';
cfgS.surfdownsample      = 10;
cfgS.funcolorlim         = 'maxabs';
cfgS.interactive         = 'yes';
ft_sourceplot(cfgS, sourceDiffIntAVG);

%% error occurs with 
% addpath('/home/common/matlab/fieldtrip/');
%% but not with 
% addpath('/home/common/matlab/fieldtrip-20131231/');

%% Error description
%
% the input is volume data with dimensions [91 109 91]
% scaling anatomy to [0 1]
% no masking parameter
% The source data is defined on a 3D grid, interpolation to a surface mesh will be performed
% Warning: could not determine filetype of surface_l4_both.mat 
% > In fileio/private/warning_once at 158
%   In ft_filetype at 1189
%   In ft_read_headshape at 110
%   In ft_sourceplot at 890 
% Error using ft_read_headshape (line 888)
% unknown fileformat "unknown" for head shape information
% 
% Error in ft_sourceplot (line 890)
%     surf  = ft_read_headshape(cfg.surffile);
