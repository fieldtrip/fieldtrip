function test_pull393

% WALLTIME 00:10:00
% MEM 2gb

%% construct some data
% data_eeg
data_eeg.label = {'eeg 1';'eeg 2';'eeg 3'};
data_eeg.trial{1,1} = randn(3,10);
data_eeg.time{1,1} = 1:10;
data_eeg.elec.label = data_eeg.label;
data_eeg.elec.elecpos = [1 1 1; 2 2 2; 3 3 3];
data_eeg.elec.chanpos = [1 1 1; 2 2 2; 3 3 3];
data_eeg.elec.tra = eye(3);

% data_eeg2 (a different EEG recording: same labels but different elecpos: this should break)
data_eeg2 = data_eeg;
data_eeg2.elec.elecpos = data_eeg.elec.elecpos+randn(1); % note the random offset
data_eeg2.elec.chanpos = data_eeg2.elec.elecpos;

% data_meg
data_meg.label = {'meg 1';'meg 2';'meg 3'};
data_meg.trial{1,1} = randn(3,10);
data_meg.time{1,1} = 1:10;
data_meg.grad.label = data_meg.label;
data_meg.grad.coilpos = [1 1 1; 2 2 2; 3 3 3];
data_meg.grad.coilori = NaN(3,3);
data_meg.grad.chanpos = [1 1 1; 2 2 2; 3 3 3];
data_meg.grad.chanori = NaN(3,3);
data_meg.grad.tra = eye(3);

% data_nirs
data_nirs.label = {'nirs 1';'nirs 2';'nirs 3'};
data_nirs.trial{1,1} = randn(3,10);
data_nirs.time{1,1} = 1:10;
data_nirs.opto.label = data_nirs.label;
data_nirs.opto.optopos = [1 1 1; 2 2 2; 3 3 3];
data_nirs.opto.chanpos = [1 1 1; 2 2 2; 3 3 3];
data_nirs.opto.tra = eye(3); 

% data_ieeg monopolar
data_ieeg.label = {'ieeg 1';'ieeg 2';'ieeg 3'};
data_ieeg.trial{1,1} = randn(3,10);
data_ieeg.time{1,1} = 1:10;
data_ieeg.elec.label = data_ieeg.label;
data_ieeg.elec.elecpos = [1 1 1; 2 2 2; 3 3 3];
data_ieeg.elec.chanpos = [1 1 1; 2 2 2; 3 3 3];
data_ieeg.elec.tra = eye(3);

% data_ieeg bipolar
data_ieegb.label = {'ieeg 1 - ieeg 2';'ieeg 2 - ieeg 3'};
data_ieegb.trial{1,1} = randn(2,10);
data_ieegb.time{1,1} = 1:10;
data_ieegb.elec.label = data_ieegb.label;
data_ieegb.elec.elecpos = [1 1 1; 2 2 2; 3 3 3];
data_ieegb.elec.chanpos = [1.5 1.5 1.5; 2.5 2.5 2.5];
data_ieegb.elec.tra = [+1 -1  0; 0 +1 -1];
data_ieegb.elec.labelold = {'ieeg 1';'ieeg 2';'ieeg 3'}; % note the labelold field due to re-referencing
data_ieegb.elec.chanposold = [1 1 1; 2 2 2; 3 3 3];                     % note the chanposold field due to re-referencing

%% do a sanity check on the bipolar iEEG electrode definition

bipolar.labelold = {
  'ieeg 1'
  'ieeg 2'
  'ieeg 3'
  };
bipolar.labelnew = {
  'ieeg 1 - ieeg 2'
  'ieeg 2 - ieeg 3'
  };
bipolar.tra = [
  1 -1 0
  0 1 -1
  ];

elecb = ft_apply_montage(data_ieeg.elec, bipolar);

assert(isequal(elecb, data_ieegb.elec));

%% convert to timlock representation

cfg = [];
timelock_eeg   = ft_timelockanalysis(cfg, data_eeg);
timelock_meg   = ft_timelockanalysis(cfg, data_meg);
timelock_nirs  = ft_timelockanalysis(cfg, data_nirs);
timelock_ieeg  = ft_timelockanalysis(cfg, data_ieeg);
timelock_ieegb = ft_timelockanalysis(cfg, data_ieegb);


%% convert to freq representation

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
freq_eeg   = ft_freqanalysis(cfg, data_eeg);
freq_meg   = ft_freqanalysis(cfg, data_meg);
freq_nirs  = ft_freqanalysis(cfg, data_nirs);
freq_ieeg  = ft_freqanalysis(cfg, data_ieeg);
freq_ieegb = ft_freqanalysis(cfg, data_ieegb);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% append raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % eeg and eeg2: this should break because of inconsistent label and chanpos
% cfg = [];
% cfg.appendsens = 'yes';
% append0 = ft_appenddata(cfg, data_eeg, data_eeg2);
% assert(isequal(numel(append0.label),3)) % 3 labels (concat across rpt dim)
% assert(isequal(numel(append0.elec.label),3)) % 3 elec labels

% eeg and meg
cfg = [];
cfg.appendsens = 'yes';
append1 = ft_appenddata(cfg, data_eeg, data_meg);
assert(isequal(numel(append1.label),6)) % 6 labels
assert(isequal(numel(append1.elec.label),3)) % 3 elec labels
assert(isequal(numel(append1.grad.label),3)) % 3 grad labels
assert(isequal(size(append1.elec.tra,1),3)) % 3 tra rows
assert(isequal(size(append1.grad.tra,1),3)) % 3 tra rows

% eeg and nirs
cfg = [];
cfg.appendsens = 'yes';
append2 = ft_appenddata(cfg, data_eeg, data_nirs);
assert(isequal(numel(append2.label),6)) % 6 labels
assert(isequal(numel(append2.elec.label),3)) % 3 elec labels
assert(isequal(numel(append2.opto.label),3)) % 3 opto labels
assert(isequal(size(append2.elec.tra,1),3)) % 3 tra rows
assert(isequal(size(append2.opto.tra,1),3)) % 3 tra rows

% monopolar and bipolar
cfg = [];
cfg.appendsens = 'yes';
append3 = ft_appenddata(cfg, data_ieeg, data_ieegb);
assert(isequal(numel(append3.label),5)) % 5 ieeg labels
assert(isequal(numel(append3.elec.label),5)) % 5 ieeg labels
assert(isequal(size(append3.elec.elecpos,1),3)) % 3 original elecpos
assert(isequal(size(append3.elec.tra,1),5)) % 5 tra rows

% do not append (discard inconsistent sens information)
append4 = ft_appenddata([], data_eeg, data_meg);
assert(isequal(numel(append4.label),6)) % 6 labels
assert(~isfield(append4, 'elec')) % no elec struc
assert(~isfield(append4, 'grad')) % no grad struc

% do not append (discard inconsistent sens information)
append5 = ft_appenddata([], data_eeg, data_nirs);
assert(isequal(numel(append5.label),6)) % 6 labels
assert(~isfield(append5, 'elec')) % no elec struc
assert(~isfield(append5, 'opto')) % no opto struc

% monopolar and bipolar
append6 = ft_appenddata([], data_ieeg, data_ieegb);
assert(isequal(numel(append6.label),5)) % 5 ieeg labels
assert(~isfield(append6, 'elec')) % no elec struc

% use the same sensor information
append7 = ft_appenddata([], data_eeg, data_eeg);
assert(isequal(numel(append7.label),3)) % 3 labels
assert(isfield(append7, 'elec')) % elec struc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% append timelock data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eeg and meg
cfg = [];
cfg.appendsens = 'yes';
append1 = ft_appendtimelock(cfg, timelock_eeg, timelock_meg);
assert(isequal(numel(append1.label),6)) % 6 labels
assert(isequal(numel(append1.elec.label),3)) % 3 elec labels
assert(isequal(numel(append1.grad.label),3)) % 3 grad labels
assert(isequal(size(append1.elec.tra,1),3)) % 3 tra rows
assert(isequal(size(append1.grad.tra,1),3)) % 3 tra rows

% eeg and nirs
cfg = [];
cfg.appendsens = 'yes';
append2 = ft_appendtimelock(cfg, timelock_eeg, timelock_nirs);
assert(isequal(numel(append2.label),6)) % 6 labels
assert(isequal(numel(append2.elec.label),3)) % 3 elec labels
assert(isequal(numel(append2.opto.label),3)) % 3 opto labels
assert(isequal(size(append2.elec.tra,1),3)) % 3 tra rows
assert(isequal(size(append2.opto.tra,1),3)) % 3 tra rows

% monopolar and bipolar
cfg = [];
cfg.appendsens = 'yes';
append3 = ft_appendtimelock(cfg, timelock_ieeg, timelock_ieegb);
assert(isequal(numel(append3.label),5)) % 5 ieeg labels
assert(isequal(numel(append3.elec.label),5)) % 5 ieeg labels
assert(isequal(size(append3.elec.elecpos,1),3)) % 3 original elecpos
assert(isequal(size(append3.elec.tra),[5 3]))

% do not append (discard inconsistent sens information)
append4 = ft_appendtimelock([], timelock_eeg, timelock_meg);
assert(isequal(numel(append4.label),6)) % 6 labels
assert(~isfield(append4, 'elec')) % no elec struc
assert(~isfield(append4, 'grad')) % no grad struc

% do not append (discard inconsistent sens information)
append5 = ft_appendtimelock([], timelock_eeg, timelock_nirs);
assert(isequal(numel(append5.label),6)) % 6 labels
assert(~isfield(append5, 'elec')) % no elec struc
assert(~isfield(append5, 'opto')) % no opto struc

% monopolar and bipolar
append6 = ft_appendtimelock([], timelock_ieeg, timelock_ieegb);
assert(isequal(numel(append6.label),5)) % 5 ieeg labels
assert(~isfield(append6, 'elec')) % no elec struc

% use the same sensor information
append7 = ft_appendtimelock([], timelock_eeg, timelock_eeg);
assert(isequal(numel(append7.label),3)) % 3 labels
assert(isfield(append7, 'elec')) % elec struc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% append freq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eeg and meg
cfg = [];
cfg.appendsens = 'yes';
append1 = ft_appendfreq(cfg, freq_eeg, freq_meg);
assert(isequal(numel(append1.label),6)) % 6 labels
assert(isequal(numel(append1.elec.label),3)) % 3 elec labels
assert(isequal(numel(append1.grad.label),3)) % 3 grad labels
assert(isequal(size(append1.elec.tra,1),3)) % 3 tra rows
assert(isequal(size(append1.grad.tra,1),3)) % 3 tra rows

% eeg and nirs
cfg = [];
cfg.appendsens = 'yes';
append2 = ft_appendfreq(cfg, freq_eeg, freq_nirs);
assert(isequal(numel(append2.label),6)) % 6 labels
assert(isequal(numel(append2.elec.label),3)) % 3 elec labels
assert(isequal(numel(append2.opto.label),3)) % 3 opto labels
assert(isequal(size(append2.elec.tra,1),3)) % 3 tra rows
assert(isequal(size(append2.opto.tra,1),3)) % 3 tra rows

% monopolar and bipolar
cfg = [];
cfg.appendsens = 'yes';
append3 = ft_appendfreq(cfg, freq_ieeg, freq_ieegb);
assert(isequal(numel(append3.label),5)) % 5 ieeg labels
assert(isequal(numel(append3.elec.label),5)) % 5 ieeg labels
assert(isequal(size(append3.elec.elecpos,1),3)) % 3 original elecpos
assert(isequal(size(append3.elec.tra,1),5)) % 5 tra rows

% do not append (discard inconsistent sens information)
append4 = ft_appendfreq([], freq_eeg, freq_meg);
assert(isequal(numel(append4.label),6)) % 6 labels
assert(~isfield(append4, 'elec')) % no elec struc
assert(~isfield(append4, 'grad')) % no grad struc

% do not append (discard inconsistent sens information)
append5 = ft_appendfreq([], freq_eeg, freq_nirs);
assert(isequal(numel(append5.label),6)) % 6 labels
assert(~isfield(append5, 'elec')) % no elec struc
assert(~isfield(append5, 'opto')) % no opto struc

% monopolar and bipolar
append6 = ft_appendfreq([], freq_ieeg, freq_ieegb);
assert(isequal(numel(append6.label),5)) % 5 ieeg labels
assert(~isfield(append6, 'elec')) % no elec struc

% use the same sensor information
append7 = ft_appendfreq([], freq_eeg, freq_eeg);
assert(isequal(numel(append7.label),3)) % 3 labels
assert(isfield(append7, 'elec')) % elec struc