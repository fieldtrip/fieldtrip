function test_shuffle_coilorder

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_prepare_headmodel ft_compute_leadfield
% DATA no

% With the commit https://github.com/fieldtrip/fieldtrip/commit/61930a70a8db47da4b5778bb20ff505d122c11c1
% the order of the coils changed in coilpos, coilori and tra. It used to be first all
% bottom coils, then all top coils, but now it is bottom-top-bottom-top-... This led
% to the concern that the OpenMEEG implementaiton might be coil-order specific,
% which is tested with this script.

%%

% create headshape points on a sphere
n = 10*1000;
pnt = randn(n,3);
headshape.pos = diag(1./sqrt(sum(pnt.^2,2))) * pnt * 0.08;
headshape.unit = 'm';

cfg = [];
cfg.method = 'headshape';
cfg.headshape = headshape;
cfg.numvertices = 1000; % downsample to make the vertices more evenly distributed
mesh = ft_prepare_mesh(cfg);

ft_plot_mesh(mesh)

%%

cfg = [];
cfg.method = 'singlesphere';
headmodel_sphere = ft_prepare_headmodel(cfg, mesh);

cfg = [];
cfg.method = 'openmeeg';
headmodel_openmeeg = ft_prepare_headmodel(cfg, mesh);

cfg = [];
cfg.method = 'singleshell';
headmodel_singleshell = ft_prepare_headmodel(cfg, mesh);

%%

grad = ft_read_sens('ctf151.mat');
dippos = [0 0 0.07];

% shuffle the coilspwd
[nchan, ncoil] = size(grad.tra);
shuffle = randperm(ncoil);

grad.coilpos = grad.coilpos(shuffle,:);
grad.coilori = grad.coilori(shuffle,:);
grad.tra = grad.tra(:,shuffle);

lf_sphere = ft_compute_leadfield(dippos, grad, headmodel_sphere);
lf_openmeeg = ft_compute_leadfield(dippos, grad, headmodel_openmeeg);

% singleshell requires the forwpar to be computed in FT_PREPARE_VOL_SENS
[headmodel_singleshell, grad] = ft_prepare_vol_sens(headmodel_singleshell, grad);
lf_singleshell = ft_compute_leadfield(dippos, grad, headmodel_singleshell);

%%

figure; hold on

dipmom = [1 0 0]';
plot(lf_sphere * dipmom, lf_openmeeg * dipmom, 'r.')
plot(lf_sphere * dipmom, lf_singleshell * dipmom, 'ro')

assert(corr(lf_sphere * dipmom, lf_openmeeg * dipmom)>0.99);
assert(corr(lf_sphere * dipmom, lf_singleshell * dipmom)>0.99);

dipmom = [0 1 0]';
plot(lf_sphere * dipmom, lf_openmeeg * dipmom, 'g.')
plot(lf_sphere * dipmom, lf_singleshell * dipmom, 'go')

dipmom = [0 0 1]';
plot(lf_sphere * dipmom, lf_openmeeg * dipmom, 'b.') % should be close to zero
plot(lf_sphere * dipmom, lf_singleshell * dipmom, 'bo')
