function test_bug3119

% WALLTIME 00:20:00
% MEM 2gb

% TEST ft_dipolefitting dipole_fit

%% load template mri
ftdir = fileparts(which('ft_defaults'));
mri = ft_read_mri(fullfile(ftdir,'template','headmodel','standard_mri.mat'));

%% used later to plot dipoles
cfg = [];
cfg.resolution = 1;
cfg.xrange = [-100 100];
cfg.yrange = [-140 140];
cfg.zrange = [-80 120];
mris = ft_volumereslice(cfg, mri);
mris = ft_convert_units(mris, 'cm');

%% load template segmented mri and rearrange fields
seg = ft_read_mri(fullfile(ftdir,'template','headmodel','standard_seg.mat'));

brain = zeros(size(seg.seg));
bri = find(seg.seg==3);
brain(bri)=1; clear bri;
skull = zeros(size(seg.seg));
ski = find(seg.seg==2);
skull(ski)=1; clear ski;
scalp = zeros(size(seg.seg));
sci = find(seg.seg==1);
scalp(sci)=1; clear sci;

seg.brain = brain; clear brain;
seg.scalp = scalp; clear scalp;
seg.skull = skull; clear skull;
seg = rmfield(seg,'seg');

vol = ft_read_vol(fullfile(ftdir,'template','headmodel','standard_bem.mat'));

%%
figure;
ft_plot_mesh(vol.bnd(1),'edgecolor', 'none', 'facecolor', 'r'); % scalp
ft_plot_mesh(vol.bnd(2),'edgecolor', 'none', 'facecolor', 'g'); % skull
ft_plot_mesh(vol.bnd(3),'edgecolor', 'none', 'facecolor', 'b'); % brain
alpha 0.3;
view(132, 14)


%% get custom electrodes
label = {'Fp1';'AF7';'AF3';'F1';'TP9';'F5h';'F7';'FT7';'FC5h';'PO9';'FC1';'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'I1';'P5h';'P7';'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';'AF4';'AFz';'Fz';'F2';'TP10';'F6h';'F8';'FT8';'FC6h';'PO10';'FC2';'FCz';'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'I2';'P6h';'P8';'P10';'PO8';'PO4';'O2'};
elec = ft_read_sens(fullfile(ftdir,'template','electrode','standard_1005.elc'));

[sel1,sel2] = match_str(label,elec.label);
% isequal(elec.label(sel2),label(sel1))

chanpos = elec.chanpos(sel2,:);
elecpos = elec.elecpos(sel2,:);
elec.chanpos = chanpos;
elec.elecpos = elecpos;
elec.label = label;

%%
figure;
hold on
ft_plot_sens(elec, 'style', 'ob','label','label');
ft_plot_vol(vol, 'facealpha', 0.5, 'edgecolor', 'none'); % "lighting phong" does not work with opacity
material dull;
camlight;

%%
cfg                 = [];
cfg.elec            = ft_convert_units(elec,'cm');
cfg.headmodel       = ft_convert_units(vol,'cm');
cfg.reducerank      = 3;
cfg.channel         = 'all';
cfg.grid.resolution = 1;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
grid = ft_prepare_leadfield(cfg);


%% load simulated data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3119.mat'));

figure
plot(datat1.time, datat1.avg)
figure
semilogy(svd(datat1.avg), '.') % this shows two sources in the data

%%
% the following was the original, which is not a correct specification of the cfg

cfg            = [];
cfg.channel    = 'all';
cfg.elec       = ft_convert_units(elec,'cm');
cfg.grid       = ft_convert_units(grid,'cm'); % note that this grid is incorrect, since Npos*3
cfg.headmodel  = ft_convert_units(vol,'cm');
cfg.senstype   = 'eeg';
cfg.latency    = [0.2 0.3];
cfg.reducerank = 3;
cfg.gridsearch = 'yes';
cfg.nonlinear  = 'yes';
cfg.model      = 'moving';
cfg.numdipoles = 2;
cfg.symmetry   = 'x';
cfg.dipfit.constr.sequential = 1;

try
  dipole = ft_dipolefitting(cfg, datat1);
  errorthrown = false;
catch
  errorthrown = true;
end
assert(errorthrown, 'the expected error was not thrown');

%% 
% let's try some other configurations that should work

cfg                 = [];
cfg.elec            = ft_convert_units(elec,'cm');
cfg.headmodel       = ft_convert_units(vol,'cm');
cfg.reducerank      = 3;
cfg.channel         = 'all';
cfg.grid.resolution = 1;
cfg.grid.unit       = 'cm';
cfg.symmetry        = 'x';
grid2 = ft_prepare_leadfield(cfg);

cfg            = [];
cfg.channel    = 'all';
cfg.elec       = ft_convert_units(elec,'cm');
cfg.grid       = ft_convert_units(grid2,'cm');
cfg.headmodel  = ft_convert_units(vol,'cm');
cfg.senstype   = 'eeg';
cfg.latency    = [0.2 0.3];
cfg.reducerank = 3;
cfg.gridsearch = 'yes';
cfg.nonlinear  = 'yes';
cfg.model      = 'moving';
cfg.numdipoles = 2;
cfg.symmetry   = 'x';

dipole = ft_dipolefitting(cfg, datat1);
  
%%

cfg            = [];
cfg.channel    = 'all';
cfg.elec       = ft_convert_units(elec,'cm');
cfg.grid       = ft_convert_units(grid,'cm');
cfg.headmodel  = ft_convert_units(vol,'cm');
cfg.senstype   = 'eeg';
cfg.latency    = [0.2 0.3];
cfg.reducerank = 3;
cfg.gridsearch = 'yes';
cfg.nonlinear  = 'yes';
cfg.model      = 'moving';

dipole = ft_dipolefitting(cfg, datat1);

%% 

cfg            = [];
cfg.channel    = 'all';
cfg.elec       = ft_convert_units(elec,'cm');
cfg.grid       = ft_convert_units(grid,'cm');
cfg.headmodel  = ft_convert_units(vol,'cm');
cfg.senstype   = 'eeg';
cfg.latency    = [0.2 0.3];
cfg.reducerank = 3;
cfg.gridsearch = 'yes';
cfg.nonlinear  = 'yes';
cfg.model      = 'regional';
cfg.numdipoles = 2;
cfg.symmetry   = 'x';

dipole = ft_dipolefitting(cfg, datat1);

%%

original = rmfield(dipole, 'dip');
original.avg = dipole.Vdata;

fitted = rmfield(dipole, 'dip');
fitted.avg = dipole.Vmodel;

cfg = [];
cfg.layout = 'elec1005.lay';
ft_multiplotER(cfg, original, fitted);



%%
% figure;
% ft_plot_slice(mris.anatomy, 'transform', mris.transform, 'location',  [-3.1   -4.6    -1.4], 'orientation', [0 1 0], 'resolution', 0.1);
% ft_plot_slice(mris.anatomy, 'transform', mris.transform, 'location',  [-3.1   -4.6    -1.4], 'orientation', [1 0 0], 'resolution', 0.1);
% ft_plot_slice(mris.anatomy, 'transform', mris.transform, 'location',  [-3.1   -4.6    -1.4], 'orientation', [0 0 1], 'resolution', 0.1);
% axis tight;
% axis off;
%
% for frame=1:17;
%   ft_plot_dipole(dipole.dip(frame).pos(1,:), dipole.dip(frame).mom(1:3,:),'color','r');
%   ft_plot_dipole(dipole.dip(frame).pos(2,:), dipole.dip(frame).mom(4:6,:),'color','b');
% end
