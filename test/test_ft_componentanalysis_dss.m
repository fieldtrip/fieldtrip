function test_ft_componentanalysis_dss(datadir)

% WALLTIME 00:20:00
% MEM 12gb
% DEPENDENCY ft_componentanalysis ft_artifact_zvalue denoise_avg2

if nargin==0
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/workshop/cuttingeeg2021/data-minimal');
else
  % assume the input argument to correspond to where the data can be found
  % locally
end

% load the eog
cfg              = [];
cfg.dataset      = fullfile(datadir, 'sub-02/eeg/sub-02_task-language_eeg.vhdr');
cfg.channel      = {'50' '64'};
cfg.reref        = 'yes';
cfg.implicitref  = []; % this is the default, we mention it here to be explicit
cfg.refchannel   = {'50'};
eogv             = ft_preprocessing(cfg);

% load the data
cfg             = [];
cfg.dataset     = fullfile(datadir, 'sub-02/eeg/sub-02_task-language_eeg.vhdr');
cfg.reref       = 'yes';
cfg.channel     = {'all' '-61' '-62' '-63' '-64'};
cfg.implicitref = 'M1';             % the implicit (non-recorded) reference channel is added to the data
cfg.refchannel  = {'M1', '53'};     % the average of these will be the new reference, note that '53' corresponds to the right mastoid (M2)
data            = ft_preprocessing(cfg);

% process the eog a bit more
cfg              = [];
cfg.channel      = '64';
eogv             = ft_selectdata(cfg, eogv);
eogv.label       = {'EOGv'};

cfg            = [];
cfg.length     = 10;
eogv_segmented = ft_redefinetrial(cfg, eogv);


% cfg            = [];
% cfg.method     = 'summary';
% cfg.keeptrial  = 'nan';
% eogv_segmented_clean = ft_rejectvisual(cfg, eogv_segmented);

% blockbreaks are approximately here
%510001-730000
%1225001-1515000

cfg = [];
cfg.trials = [1:102 147:245 304:numel(eogv_segmented.trial)];
eogv_segmented_clean = ft_selectdata(cfg, eogv_segmented);

% identify the peaks
cfg                                  = [];
cfg.artfctdef.zvalue.channel         = 'EOGv';
cfg.artfctdef.zvalue.cutoff          = 3;
cfg.artfctdef.zvalue.interactive     = 'no';%'yes';
cfg.artfctdef.zvalue.bpfilter        = 'yes';
cfg.artfctdef.zvalue.bpfreq          = [1 20];
cfg.artfctdef.zvalue.bpfilttype      = 'firws';
cfg.artfctdef.zvalue.rectify         = 'yes';
cfg.artfctdef.zvalue.artfctpeak      = 'yes';
cfg.artfctdef.zvalue.artfctpeakrange = [-.1 .3];
cfg = ft_artifact_zvalue(cfg, eogv_segmented_clean);

% put the artifact description into a params structure
params          = [];
params.demean   = true;
params.artifact = cfg.artfctdef.zvalue.artifact;

dataorig = data;

cfg        = [];
cfg.length = 5;
data       = ft_redefinetrial(cfg, dataorig);

% cfg        = [];
% cfg.method = 'summary';
% data       = ft_rejectvisual(cfg, data);

cfg                   = [];
cfg.method            = 'dss';
cfg.dss.denf.function = 'denoise_avg2';
cfg.dss.denf.params   = params;
cfg.numcomponent      = 5;
cfg.cellmode          = 'yes';
comp                  = ft_componentanalysis(cfg, data);

cfg.cellmode = 'no';
comp2        = ft_componentanalysis(cfg, data);

cfg           = []; 
cfg.layout    = 'dccn_customized_acticap64.mat';
cfg.component = 1:4;
figure;ft_topoplotIC(cfg, comp);
figure;ft_topoplotIC(cfg, comp2);

