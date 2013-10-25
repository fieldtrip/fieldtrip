function test_testbug1563(datainfo, version)

% MEM 1gb
% WALLTIME 00:03:05

% TEST test_bug1599
% TEST ft_sourceanalysis beamformer_lcmv

% fixedori is set correctly in beamformer_lcmv, the only problem is
% that sourceanalysis seems to not bother about this when using a
% predefined filter

if nargin<1
  datainfo = ref_datasets;
end
if nargin<2
  version = 'latest';
end

if ~isequal('PCWIN64', computer)
  warning('this bug only occurs under windows and matlab 64 bit')
end

%% Try to reproduce
% load your very favourite data
k = 10; % ctf data
inputfile = fullfile(datainfo(k).origdir,version,'timelock',datainfo(k).type,['timelock_trl_' datainfo(k).datatype '.mat']);
load(inputfile);

% make grid'n'vol
vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.unit = 'cm';
vol.type = 'singlesphere';

grid = [];
grid.resolution = 2.5;
grid.xgrid = 'auto';
grid.ygrid = 'auto';
grid.zgrid = 'auto';

% compute filter
cfg                 = [];
cfg.vol             = vol;
cfg.grid            = grid;
cfg.method          = 'lcmv';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   ='no'; 
source              = ft_sourceanalysis(cfg, timelock); 

if isunix
    template = '/home/common/matlab/fieldtrip/external/spm8/templates/T1.nii'; %template is distributed with spm
elseif ispc
    template = 'H:\common\matlab\fieldtrip\external\spm8\templates/T1.nii'; %template is distributed with spm
end

template_mri = ft_read_mri(template);
cfg            = [];
cfg.voxelcoord = 'no';
cfg.parameter  = {'avg.pow'};
cfg.interpmethod = 'spline';
cfg.coordsys   = 'mni';
source_int  = ft_sourceinterpolate(cfg, source, template_mri);
