function test_issue1395

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY opm_fil

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/opm_fil'));

%% the data is from https://github.com/tierneytim/OPM

channelsfile  = 'testData/channels.tsv';
coordsysfile  = 'testData/coordsystem.json';
datafile      = 'testData/meg.bin';
headerfile    = 'testData/meg.json';
positionsfile = 'testData/positions.tsv';

mrifile       = 'testData/T1w.nii';
cortexfile    = 'testData/testCustom.gii';

%%

hdr = ft_read_header(headerfile, 'headerformat', 'opm_fil');
dat = ft_read_data(headerfile, 'headerformat', 'opm_fil', 'dataformat', 'opm_fil');
event = ft_read_event(headerfile, 'headerformat', 'opm_fil', 'eventformat', 'opm_fil');

%%

cfg = [];
cfg.headerformat = 'opm_fil';
cfg.dataformat = 'opm_fil';
cfg.eventformat = 'opm_fil';
cfg.dataset = datafile;
data = ft_preprocessing(cfg);

%%

cfg = [];
cfg.headerformat = 'opm_fil';
cfg.dataformat = 'opm_fil';
cfg.eventformat = 'opm_fil';
cfg.dataset = datafile;
cfg.trialdef.eventtype = '?';
ft_definetrial(cfg); % let us look at the event types and values

cfg.trialdef.eventtype = 'TRIG1';
cfg.trialdef.eventvalue = 16;
cfg.trialdef.prestim = 0.2;
cfg.trialdef.poststim = 0.8;
cfg = ft_definetrial(cfg);

cfg.demean = 'yes';
cfg.baselinewindow = [-inf 0];
data = ft_preprocessing(cfg);

%%

% cfg = [];
% cfg.method = 'summary';
% data_clean = ft_rejectvisual(cfg, data);

% prevent user interaction in the test script
data_clean = data;

%%

cfg = [];
timelock = ft_timelockanalysis(cfg, data_clean);

%%
cfg = [];
cfg.layout = 'ordered';
cfg.channel = 1:17;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, timelock);

%%

grad = timelock.grad;

figure
ft_plot_sens(grad, 'label', 'label');
ft_plot_axes(grad);

% grad = ft_determine_coordsys(grad);
%
% Do you want to change the anatomical labels for the axes [y, N]? y
% What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? r
% What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? a
% What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? s
% Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), n(ot a landmark)? n

% prevent user interaction in the test script
grad.coordsys = 'ras'; 

%%

cfg = [];
cfg.grad = grad;
cfg.projection = 'orthographic';
cfg.viewpoint = 'topright';
cfg.outline = 'no';
layout = ft_prepare_layout(cfg);

figure
ft_plot_layout(layout);

%%

cfg = [];
cfg.layout = layout;
cfg.showlabels = 'yes';
ft_multiplotER(cfg, timelock);

