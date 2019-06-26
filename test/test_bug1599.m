function test_bug1599(datainfo, version)

% MEM 1500mb
% WALLTIME 00:10:00

% DEPENDENCY ft_sourceanalysis beamformer_lcmv

% fixedori is set correctly in beamformer_lcmv, the only problem is
% that sourceanalysis seems to not bother about this when using a
% predefined filter

if nargin<1
  datainfo = ref_datasets;
end
if nargin<2
  version = 'latest';
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
sourcemodel.xgrid = 'auto';
sourcemodel.ygrid = 'auto';
sourcemodel.zgrid = 'auto';

% compute filter
cfg                 = [];
cfg.headmodel       = vol;
cfg.sourcemodel            = grid;
cfg.method          = 'lcmv';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   ='no'; 
filter              = ft_sourceanalysis(cfg, timelock); 

cfg               = [];
cfg.headmodel     = vol;
cfg.sourcemodel          = grid;
cfg.method        = 'lcmv';
% cfg.rawtrial      = 'yes';
cfg.sourcemodel.filter   = filter.avg.filter;
cfg.lcmv.fixedori ='yes'; 
oriyes            = ft_sourceanalysis(cfg, timelock); 
cfg.lcmv.fixedori = 'no'; 
orino             = ft_sourceanalysis(cfg, timelock); 

% cfg is of course different, see above, so remove it
orino             = rmfield(orino, 'cfg'); 
oriyes            = rmfield(oriyes, 'cfg'); 

if ~isfield(oriyes.avg, 'ori') 
  error('ori not added when it should be');
end

if isfield(orino.avg, 'ori')
  error('ori added when it should not be');
end

fprintf('test_bug1599 is all fine\n');
 
