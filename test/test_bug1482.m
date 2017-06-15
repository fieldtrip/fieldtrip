function test_bug1482

% WALLTIME 00:10:00
% MEM1gb

%%
% the first section is a test script that Robert wrote

elec = [];
elec.label = {'1', '2', '3', '4'};
elec.unit = 'cm';
elec.elecpos = [
  1 1 1
  2 2 2
  3 3 3
  4 4 4
  ];

data = [];
data.label = {'1', '2', '3', '4', '5'}; % one more than elec
% data.label = {'1', '2', '3', '4'};
for i=1:10
  data.time{i} = (1:1000)/1000;
  data.trial{i} = randn(length(data.label), 1000);
end
data.elec = elec;

%%

tmpcfg = [];
tmpcfg.reref = 'yes';
tmpcfg.implicitref = [];
tmpcfg.refchannel = 'all';

cfg = [];
cfg.montage = ft_prepare_montage(tmpcfg, data);
cfg.updatesens = 'yes';
data_reref1 = ft_preprocessing(cfg, data);

%%

cfg = [];
cfg.reref = 'yes';
cfg.implicitref = [];
cfg.refchannel = 'all';
cfg.updatesens = 'yes';
data_reref2 = ft_preprocessing(cfg, data);

%%

assert(isalmostequal(data_reref1.trial{1}, data_reref2.trial{1}, 'reltol', 1e-6));
assert(isalmostequal(data_reref1.elec, data_reref2.elec, 'reltol', 1e-6));

%%

tmpcfg = [];
tmpcfg.reref = 'yes';
tmpcfg.implicitref = 'REF';
tmpcfg.refchannel = 'all';

cfg = [];
cfg.updatesens = 'yes';
cfg.montage = ft_prepare_montage(tmpcfg, data);
data_reref1 = ft_preprocessing(cfg, data);

%%

cfg = [];
cfg.reref = 'yes';
cfg.implicitref = 'REF';
cfg.refchannel = 'all';
cfg.updatesens = 'yes';
data_reref2 = ft_preprocessing(cfg, data);

%%

assert(isalmostequal(data_reref1.trial{1}, data_reref2.trial{1}, 'reltol', 1e-6));
assert(isalmostequal(data_reref1.elec, data_reref2.elec, 'reltol', 1e-6));

%%
% the following section is the test script that Arjen wrote

%% data with elec
data_eeg.label = {'eeg 1';'eeg 2';'eeg 3'};
data_eeg.trial{1,1} = randn(3,10);
data_eeg.time{1,1}  = 1:10;
data_eeg.elec.label = data_eeg.label;
data_eeg.elec.elecpos = [
  1 1 1
  2 2 2
  3 3 3
  ];
data_eeg.elec.chanpos = [
  1 1 1
  2 2 2
  3 3 3
  ];
data_eeg.elec.tra = eye(3);

% reref using cfg.montage
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg.label);
cfg.montage.labelold   = cfg.channel;
cfg.montage.labelnew   = strcat(cfg.channel(1:end-1),' - ',cfg.channel(2:end)); % 'eeg 1-eeg 2'
cfg.montage.tra        = [1 -1];
cfg.updatesens         = 'yes';
reref_eeg_montage = ft_preprocessing(cfg, data_eeg);

assert(isequal(numel(reref_eeg_montage.label),1)) % 1 bipolar rereferenced channel
assert(isequal(numel(reref_eeg_montage.elec.label),1)) % 1 bipolar rereferenced channel
assert(~isequal(reref_eeg_montage.elec.chanpos, reref_eeg_montage.elec.elecpos)) % adjusted chanpos

% reref using cfg.reref
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg.label);
cfg.reref              = 'yes';
cfg.refchannel         = 'all';
cfg.updatesens         = 'yes';
reref_eeg_reref = ft_preprocessing(cfg, data_eeg);

assert(isequal(numel(reref_eeg_reref.label),2))       % 2 common averaged rereferenced channels
assert(isequal(numel(reref_eeg_reref.elec.label),2))  % 2 common averaged rereferenced channels

assert(~isequal(reref_eeg_reref.elec.chanpos, reref_eeg_reref.elec.elecpos)) % original chanpos

% reref using cfg.reref and cfg.implicitref
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg.label);
cfg.reref              = 'yes';
cfg.refchannel         = 'all';
cfg.implicitref        = 'eeg 3';
cfg.updatesens         = 'yes';
reref_eeg_implref = ft_preprocessing(cfg, data_eeg);

assert(isequal(numel(reref_eeg_implref.label),3))       % 3 common averaged rereferenced channels
assert(isequal(numel(reref_eeg_implref.elec.label),3))  % as per https://github.com/fieldtrip/fieldtrip/pull/415#issuecomment-299952415 
assert(isequal(reref_eeg_implref.elec.chanpos, reref_eeg_implref.elec.elecpos)) % original chanpos

%% data without elec
data_eeg2.label = {'eeg 1';'eeg 2';'eeg 3'};
data_eeg2.trial{1,1} = randn(3,10);
data_eeg2.time{1,1} = 1:10;

% reref using cfg.montage
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg2.label);
cfg.montage.labelold   = cfg.channel;
cfg.montage.labelnew   = strcat(cfg.channel(1:end-1),' - ',cfg.channel(2:end)); % 'eeg 1-eeg 2'
cfg.montage.tra        = [1 -1];
cfg.updatesens         = 'yes';
reref_eeg_montage2 = ft_preprocessing(cfg, data_eeg2);

assert(isequal(numel(reref_eeg_montage2.label),1))  % 1 bipolar rereferenced channel
assert(~isfield(reref_eeg_montage2,'elec'))         % no elec

% reref using cfg.reref
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg2.label);
cfg.reref              = 'yes';
cfg.refchannel         = 'all';
cfg.updatesens         = 'yes';
reref_eeg_reref2 = ft_preprocessing(cfg, data_eeg2);

assert(isequal(numel(reref_eeg_reref2.label),2))  % 2 common averaged rereferenced channels
assert(~isfield(reref_eeg_reref2,'elec'))         % no elec

