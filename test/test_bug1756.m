function test_bug1756

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_headmodeltype ft_prepare_headmodel ft_headmodel_openmeeg

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
% ensure ft_headmodeltype works for each of the volume conduction models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(ft_headmodeltype(bemcp, 'bemcp'))
assert(ft_headmodeltype(dipoli, 'dipoli'))
assert(ft_headmodeltype(singleshell, 'singleshell'))
assert(ft_headmodeltype(concentricspheres, 'concentricspheres'))
assert(ft_headmodeltype(localspheres, 'localspheres'))
assert(ft_headmodeltype(singlesphere_eeg, 'singlesphere'))
assert(ft_headmodeltype(singlesphere_meg, 'singlesphere'))

assert( ft_headmodeltype(bemcp, 'bem'))
assert( ft_headmodeltype(dipoli, 'bem'))
assert(~ft_headmodeltype(singleshell, 'bem'))

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
assert(ft_headmodeltype(t_bemcp, 'bem'));
assert(ft_headmodeltype(t_dipoli, 'bem'));
assert(~ft_headmodeltype(t_singleshell, 'bem'));
assert(ft_headmodeltype(t_singleshell,'singleshell'))
% the following ones can still be detected properly
assert(ft_headmodeltype(t_concentricspheres, 'concentricspheres'))
assert(ft_headmodeltype(t_localspheres, 'localspheres'))
assert(ft_headmodeltype(t_singlesphere_eeg, 'singlesphere'))
assert(ft_headmodeltype(t_singlesphere_meg, 'singlesphere'))




