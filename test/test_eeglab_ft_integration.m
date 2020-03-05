function test_eeglab_ft_integration

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_dipolefitting ft_checkdata

% See http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2595

% ensure that fieldtrip/external/eeglab is on the path
ft_hastoolbox('eeglab', 1)

%%
clear all

% this contains "cfg" and "data" in FieldTrip format, generated for bugzilla #2595
load test_eeglab_ft_integration_20140529T090217.mat

% test the original reported problem
source1 = ft_dipolefitting(cfg, data);

% also do the nonlinear fit
cfg.component = [1 2 3];
cfg.nonlinear = 'yes';
source2 = ft_dipolefitting(cfg, data);

% use the optimalization toolbox (which is the default if available)
source3 = ft_dipolefitting(cfg, data);

% do not use the optimalization toolbox (which is the default if available)
cfg.dipfit.optimfun = 'fminsearch';
source4 = ft_dipolefitting(cfg, data);

%%
clear all

% this contains "cfg" and "data" in FieldTrip format, generated at the 2015 Aspet workshop
load test_eeglab_ft_integration_20150528T150257.mat

source5 = ft_dipolefitting(cfg, data);

%%
clear all

% this contains the EEG structure from eeglab_data.set/fdt used during the 2019 Aspet workshop
load(dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeglab/eeglab_data.mat'))

data      = eeglab2fieldtrip(EEG, 'preprocessing');
timelock  = eeglab2fieldtrip(EEG, 'timelockanalysis');

ft_checkdata(data, 'datatype', 'raw');
ft_checkdata(timelock, 'datatype', 'timelock');

%%
clear all

% this contains the EEG structure from eeglab_data_epochs_ica.set/fdt used during the 2019 Aspet workshop
load(dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeglab/eeglab_data_epochs_ica.mat'))

data      = eeglab2fieldtrip(EEG, 'preprocessing');
timelock  = eeglab2fieldtrip(EEG, 'timelockanalysis');
comp      = eeglab2fieldtrip(EEG, 'componentanalysis');
data1     = eeglab2fieldtrip(EEG, 'chanloc');
data2     = eeglab2fieldtrip(EEG, 'chanloc_withfid');

ft_checkdata(data, 'datatype', 'raw');
ft_checkdata(timelock, 'datatype', 'timelock');
ft_checkdata(comp, 'datatype', 'comp'); % only check the ICA topographies
ft_checkdata(comp, 'datatype', 'raw+comp'); % also check the ICA timecourses
ft_datatype_sens(data1.elec);
ft_datatype_sens(data2.elec);

