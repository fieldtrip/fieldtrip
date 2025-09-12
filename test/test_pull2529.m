function test_pull2529

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_denoise_ssp
% DATA private

% the first part of the code mainly tests whether the bookkeeping in
% elec.balance.current is still working

% the second part tests the computation and application of the SSP
% projector on OPM data

%%

elec = [];
elec.elecpos = [
  0 0 1
  0 0 2
  0 0 3
  ];
elec.label = {'eeg1', 'eeg2', 'eeg3'};
elec.balance.current = 'none';

data = [];
data.label = {'eeg1', 'eeg2', 'eeg3', 'trigger'};
data.time  = {(1:1000)/1000};
data.trial = {randn(4,1000)};
data.elec = elec;

cfg = [];
cfg.method = 'pca';
cfg.channel = {'eeg*'};
cfg.updatesens = 'yes';
eeg_comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.component = []; % keep all of them
cfg.updatesens = 'yes';
eeg_rejcomp = ft_rejectcomponent(cfg, eeg_comp);
assert(length(eeg_rejcomp.label)==3);

cfg = [];
cfg.component = []; % keep all of them
cfg.updatesens = 'yes';
data_rejcomp = ft_rejectcomponent(cfg, eeg_comp, data);
assert(length(data_rejcomp.label)==4);

%%

cfg = [];
cfg.method = 'pca';
cfg.channel = {'eeg*'};
cfg.updatesens = 'yes';
cfg.numcomponent = 3;
data_comp3 = ft_componentanalysis(cfg, data);
assert(length(data_comp3.label)==3);

cfg.numcomponent = 2;
data_comp2 = ft_componentanalysis(cfg, data);
assert(length(data_comp2.label)==2);

cfg.numcomponent = 1;
data_comp1 = ft_componentanalysis(cfg, data);
assert(length(data_comp1.label)==1);

%%
% do PCA on PCA on PCA
% the 2nd and 3rd should have identity unmixing matrices (except for sign flips)

cfg = [];
cfg.method = 'pca';
cfg.channel = {'eeg*'};
cfg.updatesens = 'yes';
data_pca1 = ft_componentanalysis(cfg, data);
cfg.channel = {'all'};
data_pca2 = ft_componentanalysis(cfg, data_pca1);
cfg.channel = {'all'};
data_pca3 = ft_componentanalysis(cfg, data_pca2);

%%
% do PCA, clean, and repeat 2x
% the 1st, 2nd and 3rd cleaned version should have rank 2, 1, 0, respectively

cfg = [];
cfg.method = 'pca';
cfg.channel = {'eeg*'};
cfg.updatesens = 'yes';
data_pca1 = ft_componentanalysis(cfg, data);

cfg = [];
cfg.component = 1; % remove one
data_rejcomp1 = ft_rejectcomponent(cfg, data_pca1);

cfg = [];
cfg.method = 'pca';
cfg.channel = {'eeg*'};
cfg.updatesens = 'yes';
data_pca2 = ft_componentanalysis(cfg, data_rejcomp1);

cfg = [];
cfg.component = 1; % remove one
data_rejcomp2 = ft_rejectcomponent(cfg, data_pca2);

cfg = [];
cfg.method = 'pca';
cfg.channel = {'eeg*'};
cfg.updatesens = 'yes';
data_pca3 = ft_componentanalysis(cfg, data_rejcomp2);

cfg = [];
cfg.component = 1; % remove one
data_rejcomp3 = ft_rejectcomponent(cfg, data_pca3);

%%

% load some data and use the same as reference data
cfg = [];
cfg.dataset = dccnpath('/project/3031000.02/test/pull2529/20250905_143322_sub-emptyroomwithHPI_file-triaxis_raw.fif');
cfg.hpfilter = 'yes';
cfg.hpfreq = 1;
data = ft_preprocessing(cfg);

cfg.channel = 'MEG';
refdata = ft_preprocessing(cfg);


%%

cfg = [];
cfg.updatesens = 'no';
data_ssp = ft_denoise_ssp(cfg, data, refdata);
assert(isequal(data_ssp.grad.balance.current, {}));

cfg = [];
data_ssp = ft_denoise_ssp(cfg, data, refdata);
assert(isequal(data_ssp.grad.balance.current, {'ssp'}));

% revert the SSP
cfg = [];
cfg.ssp = 'none';
data_none = ft_denoise_ssp(cfg, data_ssp);
assert(isempty(data_none.grad.balance.current));
