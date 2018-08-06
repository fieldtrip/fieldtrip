function inspect_bug2993

% WALLTIME 00:10:00
% MEM 1gb

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2993'));
elec_original = ft_read_sens('EEGMarkersLocalite.xml');

% elec_original = ft_determine_coordsys(elec_original);
elec_original.coordsys = 'lps';

% the electrodes should be rotated around the z-axis with 180 degrees to make them consistent 
% with the MNI coordinate system (x=right, y=nose), which is expected by default in ft_prepare_layout

figure
ft_plot_sens(elec_original, 'label', 'on');
ft_plot_axes(elec_original);

%%
cfg = [];
cfg.elec = elec_original;
cfg.method = 'interactive';
elec_realigned1 = ft_electroderealign(cfg);

assert(~isequal(elec_realigned1.chanpos, elec_original.chanpos));

%%
cfg = [];
cfg.elec = elec_original;
cfg.method = 'interactive';
% this requires a headshape, but it can be empty
cfg.headshape = [];
cfg.headshape.pos  = zeros(0,3);
cfg.headshape.tri  = zeros(0,3);
cfg.headshape.unit = 'mm';
elec_realigned2 = ft_sensorrealign(cfg);

assert(~isequal(elec_realigned2.chanpos, elec_original.chanpos));

%%
cfg = [];
cfg.individual.elec = elec_original;
cfg.template.headshape = [];
cfg.template.headshape.pos  = randn(10,3);
cfg.template.headshape.tri  = zeros(0,3);
cfg.template.headshape.unit = 'mm';
cfg.method = 'interactive';
outcfg = ft_interactiverealign(cfg);

% this returns the transformation, not the updated electrode positions
assert(~isequal(outcfg.m, eye(4)));

%%

cfg = [];
cfg.elec = elec_original;
figure; ft_plot_lay(ft_prepare_layout(cfg));

cfg = [];
cfg.elec = elec_original;
cfg.rotate = 180;
figure; ft_plot_lay(ft_prepare_layout(cfg));

cfg = [];
cfg.elec = elec_realigned1;
figure; ft_plot_lay(ft_prepare_layout(cfg));


