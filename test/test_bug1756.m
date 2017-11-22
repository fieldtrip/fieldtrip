function test_bug1756

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_voltype ft_prepare_headmodel ft_headmodel_openmeeg

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1756'));

% this contains three cumulative or overlapping BEM tissues
load seg3.mat
load grad
load elec

cfg = [];
cfg.tissue = 'scalp';
cfg.numvertices = 1000;
mesh1 = ft_prepare_mesh(cfg, seg3);

cfg = [];
cfg.tissue = {'scalp', 'skull', 'brain'};
cfg.numvertices = [300 600 900];
mesh3 = ft_prepare_mesh(cfg, seg3);

cfg = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG volume conduction models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.method = 'singlesphere'; % analytical single sphere model
singlesphere_eeg = ft_prepare_headmodel(cfg, mesh1);

cfg.method = 'concentricspheres'; % analytical concentric sphere model with up to 4 spheres
concentricspheres = ft_prepare_headmodel(cfg, mesh3);

% this one does not run on my ( = roboos) Apple Desktop computer, hence I skip it for the moment
% cfg.method = 'openmeeg' % boundary element method, based on the OpenMEEG software
% openmeeg = ft_prepare_headmodel(cfg, mesh3);

cfg.method = 'bemcp'; % boundary element method, based on the implementation from Christophe Phillips
bemcp = ft_prepare_headmodel(cfg, mesh3);

cfg.method = 'dipoli'; % boundary element method, based on the implementation from Thom Oostendorp
dipoli = ft_prepare_headmodel(cfg, mesh3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEG volume conduction models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.method = 'singlesphere'; % analytical single sphere model
singlesphere_meg = ft_prepare_headmodel(cfg, mesh1);

cfg.method = 'localspheres'; % local spheres model for MEG, one sphere per channel
cfg.grad = grad;
localspheres = ft_prepare_headmodel(cfg, mesh1);
cfg = rmfield(cfg, 'grad');

cfg.method = 'singleshell'; % realisically shaped single shell approximation, based on the implementation from Guido Nolte
singleshell = ft_prepare_headmodel(cfg, mesh1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure ft_voltype works for each of the volume conduction models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(ft_voltype(bemcp, 'bemcp'))
assert(ft_voltype(dipoli, 'dipoli'))
assert(ft_voltype(singleshell, 'singleshell'))
assert(ft_voltype(concentricspheres, 'concentricspheres'))
assert(ft_voltype(localspheres, 'localspheres'))
assert(ft_voltype(singlesphere_eeg, 'singlesphere'))
assert(ft_voltype(singlesphere_meg, 'singlesphere'))

assert( ft_voltype(bemcp, 'bem'))
assert( ft_voltype(dipoli, 'bem'))
assert(~ft_voltype(singleshell, 'bem'))

% remove the type and see whether the function can recover it
t_bemcp = rmfield(bemcp,'type');
t_dipoli = rmfield(dipoli,'type');
t_singleshell = rmfield(singleshell,'type');
t_concentricspheres = rmfield(concentricspheres,'type');
t_localspheres = rmfield(localspheres,'type');
t_singlesphere_eeg = rmfield(singlesphere_eeg,'type');
t_singlesphere_meg = rmfield(singlesphere_meg,'type');

% these three are unknown, as they cannot be distinguished
% this is where the original bug crashed on
assert(ft_voltype(t_bemcp, 'bem'));
assert(ft_voltype(t_dipoli, 'bem'));
assert(~ft_voltype(t_singleshell, 'bem'));
assert(ft_voltype(t_singleshell,'singleshell'))
% the following ones can still be detected properly
assert(ft_voltype(t_concentricspheres, 'concentricspheres'))
assert(ft_voltype(t_localspheres, 'localspheres'))
assert(ft_voltype(t_singlesphere_eeg, 'singlesphere'))
assert(ft_voltype(t_singlesphere_meg, 'singlesphere'))




