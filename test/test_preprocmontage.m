%% data with elec
data_eeg.label = {'eeg 1';'eeg 2';'eeg 3'};
data_eeg.trial{1,1} = randn(3,10);
data_eeg.time{1,1} = 1:10;
data_eeg.elec.label = data_eeg.label;
data_eeg.elec.elecpos = [1 1 1; 2 2 2; 3 3 3];
data_eeg.elec.chanpos = [1 1 1; 2 2 2; 3 3 3];
data_eeg.elec.tra = eye(3); 

% reref using cfg.montage
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg.label);
cfg.montage.labelold   = cfg.channel;
cfg.montage.labelnew   = strcat(cfg.channel(1:end-1),'-',cfg.channel(2:end)); % 'eeg 1-eeg 2'
cfg.montage.tra        = [1 -1];
reref_eeg_montage = ft_preprocessing(cfg, data_eeg);
assert(isequal(numel(reref_eeg_montage.label),1)) % 1 bipolar rereferenced channel
assert(isequal(numel(reref_eeg_montage.elec.label),1)) % 1 bipolar rereferenced channel
assert(~isequal(reref_eeg_montage.elec.chanpos, reref_eeg_montage.elec.elecpos)) % adjusted chanpos

% reref using cfg.reref
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg.label);
cfg.reref              = 'yes';
cfg.refchannel         = 'all';
reref_eeg_reref = ft_preprocessing(cfg, data_eeg);
assert(isequal(numel(reref_eeg_reref.label),2)) % 2 common averaged rereferenced channels
assert(isequal(numel(reref_eeg_reref.elec.label),2)) % 2 common averaged rereferenced channels
assert(isequal(reref_eeg_reref.elec.chanpos, reref_eeg_reref.elec.elecpos)) % original chanpos

%% data without elec
data_eeg2.label = {'eeg 1';'eeg 2';'eeg 3'};
data_eeg2.trial{1,1} = randn(3,10);
data_eeg2.time{1,1} = 1:10;

% reref using cfg.montage
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg2.label);
cfg.montage.labelold   = cfg.channel;
cfg.montage.labelnew   = strcat(cfg.channel(1:end-1),'-',cfg.channel(2:end)); % 'eeg 1-eeg 2'
cfg.montage.tra        = [1 -1];
reref_eeg_montage2 = ft_preprocessing(cfg, data_eeg2);
assert(isequal(numel(reref_eeg_montage2.label),1)) % 1 bipolar rereferenced channel
assert(~isfield(reref_eeg_montage2,'elec')) % no elec

% reref using cfg.reref
cfg                    = [];
cfg.channel            = ft_channelselection({'eeg 1','eeg 2'}, data_eeg2.label);
cfg.reref              = 'yes';
cfg.refchannel         = 'all';
reref_eeg_reref2 = ft_preprocessing(cfg, data_eeg2);
assert(isequal(numel(reref_eeg_reref2.label),2)) % 2 common averaged rereferenced channels
assert(~isfield(reref_eeg_reref2,'elec')) % no elec
