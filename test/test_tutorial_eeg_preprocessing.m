function test_tutorial_eeg_preprocessing

% WALLTIME 00:45:00
% MEM 2gb

% TEST ft_definetrial ft_preprocessing

%% defining trials
dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/preprocessing_erp/s04.vhdr');

cfg = [];
cfg.trialfun     = 'trialfun_affcog';
cfg.headerfile   = dataset;
cfg = ft_definetrial(cfg);
trl = cfg.trl;

%% pre-processing and re-referencing
cfg = [];
cfg.datafile        = dataset;
cfg.trl             = trl;

% Baseline-correction options
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.2 0];

% Fitering options
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 100;

% Re-referencing options - see explanation below
cfg.reref         = 'yes';
cfg.implicitref   = 'REF';
cfg.refchannel    = {'RM' 'REF'};

data = ft_preprocessing(cfg);

cfg = [];  % use only default options                 
ft_databrowser(cfg, data);

%% extracting the EOG signals
% EOGV channel
cfg              = [];
cfg.channel      = {'53' 'LEOG'};
cfg.reref        = 'yes';
cfg.implicitref  = [];
cfg.refchannel   = {'53'};
eogv             = ft_preprocessing(cfg, data);

% only keep one channel, and rename to eogv
cfg              = [];
cfg.channel      = 'LEOG';
eogv             = ft_selectdata(cfg, eogv); 
eogv.label       = {'eogv'};

% EOGH channel
cfg              = [];
cfg.channel      = {'57' '25'};
cfg.reref        = 'yes';
cfg.implicitref  = [];
cfg.refchannel   = {'57'};
eogh             = ft_preprocessing(cfg, data);

% only keep one channel, and rename to eogh
cfg              = [];
cfg.channel      = '25';
eogh             = ft_selectdata(cfg, eogh); 
eogh.label       = {'eogh'};


% only keep all non-EOG channels
cfg         = [];
cfg.channel = setdiff(1:60, [53, 57, 25]);              % you can use either strings or numbers as selection
data        = ft_selectdata(cfg, data); 

% append the bipolar EOGH and EOGV channel to the other selected EEG channels 
cfg = [];
data = ft_appenddata(cfg, data, eogv, eogh);

disp(data.label')


%% channel layout

cfg = [];
cfg.channel  = data.label;
cfg.layout   = 'mpi_customized_acticap64';
cfg.feedback = 'yes';
lay = ft_prepare_layout(cfg);
disp(lay)

%% Artifacts
if false
  % do not run the interactive section

  cfg        = [];
  cfg.method = 'channel';
  ft_rejectvisual(cfg, data)

  cfg = [];
  cfg.method   = 'summary';
  cfg.layout   = lay;       % this allows for plotting
  cfg.channels = [1:60];    % do not show EOG channels
  data_clean   = ft_rejectvisual(cfg, data);
else
  % just copy the data over
  data_clean = data;
end

% In this exercise we suggest that you remove 8 trials with the highest 
% variance (trial numbers 22, 42, 89, 90, 92, 126, 136 and 150).
data_clean.time([22 42 89 90 92 126 136 150])       = [];
data_clean.trial([22 42 89 90 92 126 136 150])      = [];
data_clean.sampleinfo([22 42 89 90 92 126 136 150], :) = [];
data_clean.trialinfo([22 42 89 90 92 126 136 150], :)  = [];


cfg          = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data_clean);


%% Computing an plotting the ERPs

disp(data.trialinfo')

% use ft_timelockanalysis to compute the ERPs 
cfg = [];
cfg.trials = find(data_clean.trialinfo==1);
task1 = ft_timelockanalysis(cfg, data_clean);

cfg = [];
cfg.trials = find(data_clean.trialinfo==2);
task2 = ft_timelockanalysis(cfg, data_clean);

cfg = [];
cfg.layout = lay;
cfg.interactive = 'yes';
figure;
ft_multiplotER(cfg, task1, task2)

difference     = task1;                   % copy one of the structures
difference.avg = task1.avg - task2.avg;   % compute the difference ERP

cfg = [];
cfg.layout      = lay;
cfg.interactive = 'yes';
figure;
ft_multiplotER(cfg, difference)

end

