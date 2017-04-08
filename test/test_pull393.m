%% data
% EEG
eeg.label = {'EEG 1';'EEG 2';'EEG 3'};
eeg.trial{1,1} = randn(3,10);
eeg.time{1,1} = 1:10;
eeg.elec.label = eeg.label;
eeg.elec.elecpos = [1 1 1; 2 2 2; 3 3 3];
eeg.elec.chanpos = [1 1 1; 2 2 2; 3 3 3];
eeg.elec.tra = eye(3);

% MEG
meg.label = {'MEG 1';'MEG 2';'MEG 3'};
meg.trial{1,1} = randn(3,10);
meg.time{1,1} = 1:10;
meg.grad.label = meg.label;
meg.grad.coilpos = [1 1 1; 2 2 2; 3 3 3];
meg.grad.coilori = NaN(3,3);
meg.grad.chanpos = [1 1 1; 2 2 2; 3 3 3];
meg.grad.chanori = NaN(3,3);
meg.grad.tra = eye(3);

% NIRS
nirs.label = {'NIRS 1';'NIRS 2';'NIRS 3'};
nirs.trial{1,1} = randn(3,10);
nirs.time{1,1} = 1:10;
nirs.opto.label = nirs.label;
nirs.opto.optopos = [1 1 1; 2 2 2; 3 3 3];
nirs.opto.chanpos = [1 1 1; 2 2 2; 3 3 3];
nirs.opto.tra = eye(3); 

% iEEG monopolar
ieeg.label = {'iEEG 1';'iEEG 2';'iEEG 3'};
ieeg.trial{1,1} = randn(3,10);
ieeg.time{1,1} = 1:10;
ieeg.elec.label = ieeg.label;
ieeg.elec.elecpos = [1 1 1; 2 2 2; 3 3 3];
ieeg.elec.chanpos = [1 1 1; 2 2 2; 3 3 3];
ieeg.elec.tra = eye(3);

% iEEG bipolar
ieegb.label = {'iEEG 1 - iEEG 2';'iEEG 2 - iEEG 3'};
ieegb.trial{1,1} = randn(2,10);
ieegb.time{1,1} = 1:10;
ieegb.elec.label = ieegb.label;
ieegb.elec.elecpos = [1.5 1.5 1.5; 2.5 2.5 2.5];
ieegb.elec.chanpos = [1.5 1.5 1.5; 2.5 2.5 2.5];
ieegb.elec.tra = eye(2);
ieegb.elec.labelold = {'iEEG 1';'iEEG 2';'iEEG 3'}; % note the labelold field due to re-referencing
ieegb.elec.chanposold = [1 1 1; 2 2 2; 3 3 3]; % note the chanposold field due to re-referencing

%% append
% EEG and MEG
cfg = [];
cfg.appendsens = 'yes';
data1 = ft_appenddata(cfg, eeg, meg);
assert(isequal(numel(data1.label),6)) % 6 labels
assert(isequal(numel(data1.elec.label),3)) % 3 elec labels
assert(isequal(numel(data1.grad.label),3)) % 3 grad labels

% EEG and NIRS
cfg = [];
cfg.appendsens = 'yes';
data2 = ft_appenddata(cfg, eeg, nirs);
assert(isequal(numel(data2.label),6)) % 6 labels
assert(isequal(numel(data2.elec.label),3)) % 3 elec labels
assert(isequal(numel(data2.opto.label),3)) % 3 opto labels

% monopolar iEEG and bipolar iEEG
cfg = [];
cfg.appendsens = 'yes';
data3 = ft_appenddata(cfg, ieeg, ieegb);
assert(isequal(numel(data3.label),5)) % 5 ieeg labels
assert(isequal(numel(data3.elec.label),5)) % 5 ieeg labels
assert(isequal(numel(data3.elec.labelold),6)) % 6 ieeg labelolds

%% do not append (discard inconsistent sens information) 
% EEG and MEG
data4 = ft_appenddata([], eeg, meg);
assert(isequal(numel(data4.label),6)) % 6 labels
assert(~isfield(data4, 'elec')) % no elec struc
assert(~isfield(data4, 'grad')) % no grad struc

% EEG and NIRS
data5 = ft_appenddata([], eeg, nirs);
assert(isequal(numel(data5.label),6)) % 6 labels
assert(~isfield(data5, 'elec')) % no elec struc
assert(~isfield(data5, 'opto')) % no opto struc

% monopolar iEEG and bipolar iEEG
data6 = ft_appenddata([], ieeg, ieegb);
assert(isequal(numel(data6.label),5)) % 5 ieeg labels
assert(~isfield(data6, 'elec')) % no elec struc
