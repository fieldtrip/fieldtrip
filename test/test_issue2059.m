function test_issue2059

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_dipolesimulation ft_compute_leadfield

%%
% randomply distribute some electrodes overt the upper hemisphere

pos = randn(64,3);
pos(:,3) = abs(pos(:,3));
pos = diag(1./sqrt(sum(pos.^2,2))) * pos;

elec_m  = [];
elec_m.elecpos = pos*0.1;
elec_m.label = arrayfun(@num2str, 1:64, 'UniformOutput', false)';
elec_m.unit = 'm';

elec_dm  = [];
elec_dm.elecpos = pos*1;
elec_dm.label = arrayfun(@num2str, 1:64, 'UniformOutput', false)';
elec_dm.unit = 'dm';

elec_cm  = [];
elec_cm.elecpos = pos*10;
elec_cm.label = arrayfun(@num2str, 1:64, 'UniformOutput', false)';
elec_cm.unit = 'cm';

elec_mm  = [];
elec_mm.elecpos = pos*100;
elec_mm.label = arrayfun(@num2str, 1:64, 'UniformOutput', false)';
elec_mm.unit = 'mm';

%%
% construct different headmodels
% human CSF at 37 °C has a conductivity of 1.79 S/m, see https://doi.org/10.1007/978-981-10-9035-6_142

headmodel_m = [];
headmodel_m.o = [0 0 0];
headmodel_m.r = 0.10;
headmodel_m.cond = 0.33; % S/m
headmodel_m.unit = 'm';

headmodel_dm = [];
headmodel_dm.o = [0 0 0];
headmodel_dm.r = 1;
headmodel_dm.cond = 0.33/10; % S/dm
headmodel_dm.unit = 'dm';

headmodel_cm = [];
headmodel_cm.o = [0 0 0];
headmodel_cm.r = 10;
headmodel_cm.cond = 0.33/100; % S/cm
headmodel_cm.unit = 'cm';

headmodel_mm = [];
headmodel_mm.o = [0 0 0];
headmodel_mm.r = 100;
headmodel_mm.cond = 0.33/1000; % S/mm
headmodel_mm.unit = 'mm';

%%

sourcemodel_m = [];
sourcemodel_m.pos = [0 0 0.08];
sourcemodel_m.unit = 'm';

sourcemodel_dm = [];
sourcemodel_dm.pos = [0 0 0.8];
sourcemodel_dm.unit = 'dm';

sourcemodel_cm = [];
sourcemodel_cm.pos = [0 0 8];
sourcemodel_cm.unit = 'cm';

sourcemodel_mm = [];
sourcemodel_mm.pos = [0 0 80];
sourcemodel_mm.unit = 'mm';

%%
% compute the potential in Volt for a 100 nAm dipole

dipole_m =  [1 0 0]' * 100e-9; % nA*m
dipole_dm = [1 0 0]' * 1000e-9; % nA*dm
dipole_cm = [1 0 0]' * 10000e-9; % nA*cm
dipole_mm = [1 0 0]' * 100000e-9; % nA*mm

lf_m  = ft_compute_leadfield(sourcemodel_m.pos, elec_m, headmodel_m)    * dipole_m;
lf_dm = ft_compute_leadfield(sourcemodel_dm.pos, elec_dm, headmodel_dm) * dipole_dm;
lf_cm = ft_compute_leadfield(sourcemodel_cm.pos, elec_cm, headmodel_cm) * dipole_cm;
lf_mm = ft_compute_leadfield(sourcemodel_mm.pos, elec_mm, headmodel_mm) * dipole_mm;

% compute the RMS value over all electrodes
rms_m  = sqrt(mean(lf_m.^2));
rms_dm = sqrt(mean(lf_dm.^2));
rms_cm = sqrt(mean(lf_cm.^2));
rms_mm = sqrt(mean(lf_mm.^2));

assert(isalmostequal(rms_m, rms_dm, 'reltol', 100*eps));
assert(isalmostequal(rms_m, rms_dm, 'reltol', 100*eps));
assert(isalmostequal(rms_m, rms_dm, 'reltol', 100*eps));


%%
% now do the same thing with ft_dipolesimulation

cfg = [];
cfg.headmodel   = headmodel_m;
cfg.elec        = elec_m;
cfg.sourcemodel = sourcemodel_m;
cfg.sourcemodel.mom = [1 0 0]';
cfg.sourcemodel.signal{1} = 100e-9; % nA*m
data_m = ft_dipolesimulation(cfg);

cfg = [];
cfg.headmodel   = headmodel_dm;
cfg.elec        = elec_dm;
cfg.sourcemodel = sourcemodel_dm;
cfg.sourcemodel.mom = [1 0 0]';
cfg.sourcemodel.signal{1} = 1000e-9; % nA*dm
data_dm = ft_dipolesimulation(cfg);

cfg = [];
cfg.headmodel   = headmodel_cm;
cfg.elec        = elec_cm;
cfg.sourcemodel = sourcemodel_cm;
cfg.sourcemodel.mom = [1 0 0]';
cfg.sourcemodel.signal{1} = 10000e-9; % nA*cm
data_cm = ft_dipolesimulation(cfg);

cfg = [];
cfg.headmodel   = headmodel_mm;
cfg.elec        = elec_mm;
cfg.sourcemodel = sourcemodel_mm;
cfg.sourcemodel.mom = [1 0 0]';
cfg.sourcemodel.signal{1} = 100000e-9; % nA*mm
data_mm = ft_dipolesimulation(cfg);

% compute the RMS value over all electrodes
rms_m  = sqrt(mean(data_m.trial{1}.^2));
rms_dm = sqrt(mean(data_dm.trial{1}.^2));
rms_cm = sqrt(mean(data_cm.trial{1}.^2));
rms_mm = sqrt(mean(data_mm.trial{1}.^2));

assert(isalmostequal(rms_m, rms_dm, 'reltol', 100*eps));
assert(isalmostequal(rms_m, rms_dm, 'reltol', 100*eps));
assert(isalmostequal(rms_m, rms_dm, 'reltol', 100*eps));

%%
% let us now mix the units

cfg = [];
cfg.headmodel   = headmodel_m;
cfg.elec        = elec_cm;
cfg.sourcemodel = sourcemodel_mm;
% cfg.sourcemodel.pos = cfg.sourcemodel.pos/2; % this moves the dipole inward
cfg.sourcemodel.mom = [1 0 0]';
cfg.sourcemodel.signal{1} = 1;
data_mix = ft_dipolesimulation(cfg);

x = data_m.trial{1};
y = data_mix.trial{1};

figure
plot(x, y, '.');
hold on;
plot(x, x, 'r-')

% the units are not consistent, but the topography should have a perfect correlation
assert(isalmostequal(corr(x, y), 1, 'reltol', 100*eps));
