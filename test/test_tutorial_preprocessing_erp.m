function test_tutorial_preprocessing_erp

% MEM 1500mb
% WALLTIME 00:10:00

% This is a test script for the EEG preprocessing example page that Irina made for the Cracow workshop.
% See http://www.fieldtriptoolbox.org/tutorial/preprocessing_erp
%
% This version corresponds to the wiki on 24 June 2012

interactive = false;  % use this for running it without user interaction 
% interactive = true;  % use this for running it with user interaction 

datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/preprocessing_erp');
cd(datadir)

cfg = [];
cfg.trialfun     = 'trialfun_affcog';
cfg.dataset      = 's04.vhdr';
cfg = ft_definetrial(cfg);

% Baseline-correction options
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.2 0];

% Fitering options
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 100;

% Re-referencing options - see explanation below
cfg.reref         = 'yes';
cfg.refchannel    = {'RM' 'REF'};
cfg.implicitref   =  'REF';

data = ft_preprocessing(cfg);

cfg         = [];
ft_databrowser(cfg, data);

cfg         = [];
cfg.dataset = 's04.vhdr';
ft_databrowser(cfg);

load mpi_customized_acticap64.mat
disp(lay)

cfg = [];
cfg.layout = 'mpi_customized_acticap64.mat'; % alternatively you can do cfg.layout = lay
ft_layoutplot(cfg, data);

% EOGV channel
cfg              = [];
cfg.channel      = {'53' 'LEOG'};
cfg.reref        = 'yes';
cfg.implicitref  = [];
cfg.refchannel   = {'53'};
eogv             = ft_preprocessing(cfg, data);

% rename the 2nd channel to eogv
eogv.label       = {'dummy2', 'eogv'};

cfg              = [];
cfg.channel      = 'eogv';
eogv             = ft_selectdata(cfg, eogv); % only keep one channel, discard the dummy

% EOGH channel
cfg              = [];
cfg.channel      = {'57' '25'};
cfg.reref        = 'yes';
cfg.implicitref  = [];
cfg.refchannel   = {'57'};
eogh             = ft_preprocessing(cfg, data);

% rename the 1st channel to eogh
eogh.label       = {'eogh' 'dummy2'};

cfg              = [];
cfg.channel      = 'eogh';
eogh             = ft_selectdata(cfg, eogh); % only keep one channel, discard the dummy

cfg         = [];
cfg.channel = 1:60;              % you can use either strings or numbers as selection
data = ft_selectdata(cfg, data); % only keep channels 1 up to 60

% append the EOGH and EOGV channel to the 60 selected EEG channels 
cfg = [];
data = ft_appenddata(cfg, data, eogv, eogh);

disp(data.label')

if interactive
% this should not be included in automated test script execution 
cfg         = [];
cfg.method  = 'channel';
ft_rejectvisual(cfg, data);
end

if interactive
% this should not be included in automated test script execution 
cfg = [];
cfg.method   = 'summary';
cfg.layout   = lay;       % allows for plotting
cfg.channels = [1:60];    % do not show EOG channels
data_clean   = ft_rejectvisual(cfg, data);

cfg          = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data_clean);

% use ft_timelockanalysis to compute the ERPs 
cfg = [];
cfg.trials = find(data_clean.trialinfo==1);
task1= ft_timelockanalysis(cfg, data_clean);

cfg = [];
cfg.trials = find(data_clean.trialinfo==2);
task2= ft_timelockanalysis(cfg, data_clean);

cfg = [];
cfg.layout = lay;
cfg.interactive = 'yes';
ft_multiplotER(cfg, task1, task2)

difference     = task1;                   % copy one of the structures
difference.avg = task1.avg - task2.avg;   % compute the difference ERP

cfg = [];
cfg.layout = lay;
cfg.interactive = 'yes';
ft_multiplotER(cfg, difference)
end

