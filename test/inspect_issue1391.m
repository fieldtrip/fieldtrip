function inspect_issue1391

% WALLTIME 00:05:00
% MEM 2gb
% DEPENDENCY ft_rejectvisual rejectvisual_summary rejectvisual_trial rejectvisual_channel

%% create data
numtrl = 52;
numchan = 26;

data = [];
for chanlop=1:numchan
  data.label{chanlop} = char(96+chanlop); % a = 97
end
data.trialinfo = (1:numtrl)';

for trlop=1:numtrl
  % ensure the channels are distinguishable by their data (helps in the GUI)
  dat = randn(numchan,1000);
  for chanlop=1:numchan
    dat(chanlop,:) = chanlop + dat(chanlop,:);
  end
  
  % data.trial = {dat.^1 dat.^2 dat.^3 dat.^4};
  data.trial{trlop} = dat;
end

data.time = repmat({1:1000}, 1, numtrl);

%% summary mode

cfg = [];
cfg.method = 'summary';
data_clean = ft_rejectvisual(cfg, data);

% check some alternative ways of removing trials and channels, using data_clean.cfg
checkalternatives(data, data_clean);

%% trial mode, showing one trial at a time

cfg = [];
cfg.method = 'trial';
data_clean = ft_rejectvisual(cfg, data);

% check some alternative ways of removing trials and channels, using data_clean.cfg
checkalternatives(data, data_clean);

%% channel mode, showing one channel at a time

cfg = [];
cfg.method = 'channel';
data_clean = ft_rejectvisual(cfg, data);

% check some alternative ways of removing trials and channels, using data_clean.cfg
checkalternatives(data, data_clean);

%%

cfg = [];
cfg.method = 'summary';
cfg.metric = 'max';

% this one was showing problems
cfg.channel = {'a' 'c' 'd'};
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.label(:), cfg.channel(:)));

% this one was showing problems
cfg.channel = [1 3 4];
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.label(:), {'a' 'c' 'd'}'));

% try some others as well
cfg.channel = {'a' 'b' 'c' 'd'};
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.label(:), cfg.channel(:)));

% try some others as well
cfg.channel = {'d' 'b'};
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.label(:), sort(cfg.channel(:))));

% try some others as well
cfg.channel = {'b' 'c' 'd'};
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.label(:), cfg.channel(:)));

% try some others as well
cfg.channel = {'a' 'b' 'c'};
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.label(:), cfg.channel(:)));

%% trials in summary mode

cfg = [];
cfg.method = 'summary';
cfg.metric = 'max';

% try some others as well
cfg.trials = [1 2 4];
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.trialinfo(:), cfg.trials(:)));

% try some others as well
cfg.trials = [1 2 3];
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.trialinfo(:), cfg.trials(:)));

% try some others as well
cfg.trials = [2 3 4];
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.trialinfo(:), cfg.trials(:)));

% try some others as well
cfg.trials = [4 3 1];
data_clean = ft_rejectvisual(cfg, data); % don't make any additional selection, just quit
assert(isequal(data_clean.trialinfo(:), sort(cfg.trials(:))));


function checkalternatives(data, data_clean)

fprintf('-------------- alternative 1+2 --------------\n');
% this is one alternative
cfg1 = [];
cfg1.artfctdef = data_clean.cfg.artfctdef;
data_alternative1 = ft_rejectartifact(cfg1, data);
cfg2 = [];
cfg2.channel = data_clean.cfg.channel;
data_alternative2 = ft_selectdata(cfg2, data_alternative1);
assert(isequal(data_clean.trial, data_alternative2.trial));

fprintf('-------------- alternative 3 --------------\n');
% this is another alternative
cfg3 = [];
cfg3.trials = data_clean.cfg.trials;
cfg3.channel = data_clean.cfg.channel;
data_alternative3 = ft_selectdata(cfg3, data);
assert(isequal(data_clean.trial, data_alternative3.trial));
