function test_ft_electroderealign_headshape

% MEM 2gb
% WALLTIME 00:20:00

load('SubjectUCI29_elec_tal_f.mat', 'elec_tal_f'); % FIXME: please redirect to correct location at FTP/DCCN cluster
load('SubjectUCI29_hull_lh.mat', 'mesh'); % FIXME: please redirect to correct location at FTP/DCCN cluster

% keepunused
elec_tal_fr = elec_tal_f;
grids = {'LPG*', 'LTG*'};
for g = 1:numel(grids)
  cfg             = [];
  cfg.channel     = grids{g};
  cfg.keepunused  = 'yes';            % <-
  cfg.elec        = elec_tal_fr;
  cfg.method      = 'headshape';
  cfg.headshape   = mesh;
  cfg.warp        = 'dykstra2012';
  elec_tal_fr = ft_electroderealign(cfg);
end
assert(isequal(numel(elec_tal_fr.label), 152)); % realigned + unused is 152 elecs

% do not keepunused
elec_tal_fr = elec_tal_f;
grids = {'LPG*'};
cfg             = [];
cfg.channel     = grids{1};
cfg.keepunused  = 'no';             % <-
cfg.elec        = elec_tal_fr;
cfg.method      = 'headshape';
cfg.headshape   = mesh;
cfg.warp        = 'dykstra2012';
elec_tal_fr = ft_electroderealign(cfg);
assert(isequal(numel(elec_tal_fr.label), 64)); % realigned LPG has 64 elecs
