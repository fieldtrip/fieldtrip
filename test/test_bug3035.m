function test_bug3035

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_rejectcomponent ft_apply_montage
% DATA private

%% load the data

load(dccnpath('/project/3031000.02/test/bug3035/bug.mat'));

% data consists of nan-free EEG data, plus some eye tracker channels that
% contain nans. The comp structure was generated based only on the nan-free
% EEG data.

% set error function; change to warning to avoid breaking on errors
%errfun = @warning;
errfun = @error;

%% preliminary checks on data

% which channels contain nans in original data?
tmp = cat(2, data.trial{:});
hasnan = any(isnan(tmp), 2);
fprintf('the following channels have nans in the original data:\n');
disp(data.label(hasnan));

% are any of those present in comp.topolabel?
incomp = numel(match_str(data.label(hasnan), comp.topolabel));
fprintf('%d of these is/are present in comp.topolabel\n', incomp);
if incomp ~= 0
  warning('the prerequisite for this bug is not satisfied; note that this is not a failure of the test script but means that the bug does not apply to the data provided');
  return;
end

%% rejectcomponent with 3rd input argument

% NOTE: this is the original bug

cfg = [];
cfg.component = 1;
data_reject1 = ft_rejectcomponent(cfg, comp, data);

% the nans should not spread beyond the eye channels
tmp = data_reject1.trial{1}(1:42,:);
if any(isnan(tmp(:)))
  errfun('cleaned data contains nans')
end

%% rejectcomponent with 3rd input argument after zeroing nans

% NOTE: the original bug did not occur here

data1b = data;
for k = 1:numel(data1b.trial)
  data1b.trial{k}(isnan(data1b.trial{k})) = 0;
end

cfg = [];
cfg.component = 1;
data_reject1b = ft_rejectcomponent(cfg, comp, data1b);

% check
tmp = cat(1, data_reject1b.trial{:});
if any(isnan(tmp(:)))
  errfun('cleaned data contains nans');
end

%% rejectcomponent without 3rd input argument

% NOTE: the original bug did not occur here

cfg = [];
cfg.component = 1;
data_reject2 = ft_rejectcomponent(cfg, comp);

% check
tmp = cat(1, data_reject2.trial{:});
if any(isnan(tmp(:)))
  errfun('cleaned data contains nans');
end

%% rejectcomponent with 3rd input argument consisting of ONLY EEG data

% NOTE: the original bug did not occur here

% subselect only EEG data
cfg = [];
cfg.channel = 'EEG';
data_eeg = ft_selectdata(cfg, data);

cfg = [];
cfg.component = 1;
data_reject3 = ft_rejectcomponent(cfg, comp, data_eeg);

% check
tmp = cat(1, data_reject3.trial{:});
if any(isnan(tmp(:)))
  errfun('cleaned data contains nans');
end

%% if we zero-out nans beforehand, does that fix the bug?

% NOTE: the original bug did not occur here

cfg = [];
cfg.channel = 'EEG';
data_reject1c = ft_selectdata(cfg, data_reject1b);

d1 = cat(1, data_reject3.trial{:});
d2 = cat(1, data_reject1c.trial{:});

assert(isalmostequal(d1,d2,'reltol',1e-6));

