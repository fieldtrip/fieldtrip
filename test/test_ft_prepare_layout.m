function test_ft_prepare_layout

% WALLTIME 00:10:00
% MEM 4gb

% TEST ft_prepare_layout ft_plot_lay ft_plot_sens

%%
% Most of the test cases require visual inspection of the results. Nevertheless, it is useful to run them in a test batch.
% The test cases taht require manual interaction should not run in batch mode.
interactive = false;

% this corresponds to the ftp directory
cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/layout'));

%%
% make it from a figure

if interactive
  cfg = [];
  cfg.image = 'easycap_m10_equidistant61chan.png';
  cfg.bw = 'yes';
  
  layout = ft_prepare_layout(cfg);
  figure; ft_plot_lay(layout)
end

%%
% simple projection of EEG electrodes
% channel 35 is at Fpz
% the nose is along +Y

clear all
close all

elec = ft_read_sens('easycap-M10.txt');

figure
ft_plot_sens(elec, 'label', 'label')
ft_plot_axes(elec, 'fontcolor', 'm');

cfg = [];
cfg.elec = elec;

cfg.skipscale = 'yes';
cfg.skipcomnt = 'yes';
cfg.rotate = 0; % 0 is appropriate for this set of electrodes

layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

cfg.outline = 'convex';
cfg.mask    = 'convex';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

cfg.outline = 'circle';
cfg.mask    = 'convex';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

cfg.outline = 'convex';
cfg.mask    = 'circle';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

%%
% simple projection of CT151 gradiometers
% the nose is along +X

clear all
close all

load ctf151.mat % the mat files contains the variable "sens"

figure
ft_plot_sens(sens, 'label', 'label', 'chantype', 'meggrad')
ft_plot_axes(sens, 'fontcolor', 'm');

cfg = [];
cfg.grad = sens;
cfg.channel = 'MEG';

cfg.projection = 'stereographic';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.projection)

cfg.projection = 'orthographic';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.projection)

cfg.projection = 'polar';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.projection)

cfg.projection = 'gnomic'; % this works technically, but does not result in a nice layout
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.projection)

%%
% simple projection of CT151 gradiometers
% the nose is along +X

clear all
close all

load ctf151.mat % the mat files contains the variable "sens"

figure
ft_plot_sens(sens, 'label', 'label', 'chantype', 'meggrad')
ft_plot_axes(sens, 'fontcolor', 'm');

cfg = [];
cfg.grad = sens;
cfg.channel = 'MEG';

cfg.skipscale = 'yes';
cfg.skipcomnt = 'yes';

layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

cfg.outline = 'convex';
cfg.mask    = 'convex';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

cfg.outline = 'circle';
cfg.mask    = 'convex';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

cfg.outline = 'convex';
cfg.mask    = 'circle';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout)

%%
% simple projection of CT151 gradiometers
% the nose is along +X

clear all
close all

load ctf151.mat % the mat files contains the variable "sens"

figure
ft_plot_sens(sens, 'label', 'label', 'chantype', 'meggrad')
ft_plot_axes(sens, 'fontcolor', 'm');

cfg = [];
cfg.grad = sens;
cfg.channel = 'MEG';
cfg.projection = 'orthographic';

cfg.viewpoint = 'superior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint);

cfg.viewpoint = 'inferior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint);

cfg.viewpoint = 'anterior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint);

cfg.viewpoint = 'posterior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint);

cfg.channel = 'ML*';
cfg.viewpoint = 'left';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint);

cfg.channel = 'MR*';
cfg.viewpoint = 'right';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint);

%%
% WORKS if elec/grad and headshape contain the same coordsys
% WORKS for most common coordsys options

clear all
close all

load ctf151.mat % the mat files contains the variable "sens"
headshape = ft_read_headshape('Subject01.shape');
headshape.coordsys = 'ctf';

figure
ft_plot_sens(sens, 'label', 'label', 'chantype', 'meggrad')
ft_plot_axes(sens, 'fontcolor', 'm');
ft_plot_mesh(headshape)

cfg = [];
cfg.grad = sens;
cfg.channel = 'MEG';
cfg.headshape = headshape;
cfg.projection = 'orthographic';

cfg.viewpoint = 'superior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'inferior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'left';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'right';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'posterior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'anterior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)


%%
% WORKS if elec/grad and mri contain the same coordsys AND mri either has brain field or can be segmented on the fly (other segmentations can be achieved via ft_prepare_mesh and cfg.headshape)
% WORKS for most common coordsys options

clear all
close all

load ctf151.mat % the mat files contains the variable "sens"
load segmentedmri

figure
ft_plot_ortho(segmentedmri.white, 'transform', segmentedmri.transform, 'style', 'intersect', 'unit', 'mm')
ft_plot_sens(sens, 'label', 'label', 'chantype', 'meggrad', 'unit', 'mm', 'edgecolor', 'y', 'fontcolor', 'g')
ft_plot_axes([], 'unit', 'mm', 'fontcolor', 'm');

cfg = [];
cfg.grad = sens;
cfg.channel = 'MEG';
cfg.mri = segmentedmri;
cfg.projection = 'orthographic';

cfg.viewpoint = 'superior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'inferior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'posterior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.viewpoint = 'anterior';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.channel = 'ML*';
cfg.viewpoint = 'left';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.channel = 'MR*';
cfg.viewpoint = 'right';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

%%
% ok, here we actually start looking at what all the changes were about

lh = load('SubjectUCI29_cortex_lh');
rh = load('SubjectUCI29_cortex_rh');
elec = ft_read_sens('SubjectUCI29_elec_tal_f.mat');

% shift it a bit, otherwise it seems too high (!?)
elec.chanpos(:,3) = elec.chanpos(:,3) - 5;
elec.elecpos(:,3) = elec.elecpos(:,3) - 5;

% note that SubjectUCI29_elec_tal_fr is better, since projected onto the
% cortical sheet, but it lacks the depth electrodes

figure
ft_plot_sens(elec, 'label', 'label')
ft_plot_mesh(lh.mesh);
ft_plot_mesh(rh.mesh);
camlight

cfg = [];
cfg.elec = elec;
cfg.projection = 'orthographic';

cfg.headshape = lh.mesh;
cfg.channel = 'LP*';
cfg.viewpoint = 'left';
cfg.mask = 'convex';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

cfg.headshape(1) = lh.mesh;
cfg.headshape(2) = rh.mesh;
cfg.channel = {'LA*', 'LH*', 'RO*', 'RT*'};
cfg.viewpoint = 'inferior';
cfg.mask = 'none';
layout = ft_prepare_layout(cfg);
figure; ft_plot_lay(layout); title(cfg.viewpoint)

%%
% test the non-topographic layouts

clear all
close all

elec = ft_read_sens('easycap-M10.txt');

style = {'ordered', 'vertical', 'horizontal', 'butterfly', 'circular', '1column', '2column', '3column', '4column', '5column', '6column', '7column', '8column', '9column', '1row', '2row', '3row', '4row', '5row', '6row', '7row', '8row', '9row'};

for i=1:numel(style)
  cfg = [];
  cfg.channel = elec.label;
  cfg.layout = style{i};
  layout = ft_prepare_layout(cfg);
  figure; ft_plot_lay(layout); title(style{i})
end
