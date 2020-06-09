function test_bug3353

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourceparcellate test_bug3353

%% file path
datapath = dccnpath('/home/common/matlab/fieldtrip/data/test/');
filename = 'bug3353'; % contains variable source_avg generated as per:
% this is from the following tutorial http://www.fieldtriptoolbox.org/tutorial/beamformer_lcmv?s[]=beamformer&s[]=tutorial&s[]=evoked

load(fullfile(datapath, filename));

%% get atlas
ft_path       = which('ft_defaults');
[ft_path,f,e] = fileparts(ft_path);
filename      = fullfile(ft_path, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');
atlas         = ft_read_atlas(filename);

%% interpolate atlas onto source
% an aside: the documentation for ft_sourceinterpolate makes it appear as
% if only "source" can be
% interpolated onto anatomy. The function seems to be much more general

cfg              = [];
cfg.parameter    = 'tissue';
cfg.interpmethod = 'nearest';
interp_atlas     = ft_sourceinterpolate(cfg, atlas, source_avg);
interp_atlas.pos = source_avg.pos; % this is necessary due to another bug that I am also documenting later on

%% retrieve parcellations (SUCCEEDS)
cfg              = [];
cfg.method       = 'mean';
cfg.parcellation = 'tissue';
cfg.parameter    = 'pow';
parc             = ft_sourceparcellate(cfg, source_avg, interp_atlas); %% this works fine

%% retrive parcellations (used to FAIL, but has been fixed by making ft_sourceparcellate nan transparent)
cfg              = [];
cfg.method       = 'maxabs';
cfg.parcellation = 'tissue';
cfg.parameter    = 'pow';
parc             = ft_sourceparcellate(cfg, source_avg, interp_atlas); %% this fails with the following error:

%Assignment has more non-singleton rhs dimensions than non-singleton subscripts
%
%Error in ft_sourceparcellate (line 300)
%          tmp(j,:) = arraymaxabs1(dat(seg==j,:));

% it does so for max, min  and maxabs, due to max and min returning empty
% arrays when applied to empty arrays. mean and median on the other hand
% return a NaN

%% second bug (SMALL NUMERICAL INACCURACIES AREN'T "IGNORED")
% do the interpolation again, don't add pos from source_avg this time
% around
cfg              = [];
cfg.parameter    = 'tissue';
cfg.interpmethod = 'nearest';
interp_atlas     = ft_sourceinterpolate(cfg, atlas, source_avg);

%% retrieve parcellation (FAILS when no tolerance is allowed for numerical diffences)
cfg              = [];
cfg.method       = 'mean';
cfg.parcellation = 'tissue';
cfg.parameter    = 'pow';
parc             = ft_sourceparcellate(cfg, source_avg, interp_atlas); %% fails with the following error:

% Error using ft_notification (line 340)
% the source positions are not consistent with the parcellation, please use FT_SOURCEINTERPOLATE
% 
% Error in ft_error (line 39)
%   ft_notification(varargin{:});
% 
% Error in ft_sourceparcellate (line 109)
%   ft_error('the source positions are not consistent with the parcellation, please use FT_SOURCEINTERPOLATE');
 
% this is due to small numerical inaccuracies occurring at the
% interpolating level
