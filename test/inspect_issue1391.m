function inspect_issue1391

% WALLTIME 00:05:00
% MEM 2gb

% DEPENDENCY ft_rejectvisual rejectvisual_summary

%% create data

data = [];
data.label = {'a' 'b' 'c' 'd'}';
data.trialinfo = (1:4)';

% ensure the channels are distinguishable by their data (helps in the
% GUI)
dat = ones(4,1000);
dat(2,:) = 2*dat(2,:);
dat(3,:) = 3*dat(3,:);
dat(4,:) = 4*dat(4,:);

data.trial = {dat dat dat dat};
data.time = {1:1000 1:1000 1:1000 1:1000};

%% check it

cfg = [];
cfg.channel = {'a' 'c' 'd'};
cfg.trials = [1 4];
cfg.method = 'summary';
cfg.metric = 'max';
datout = ft_rejectvisual(cfg, data);

% Now I select the second channel to reject, which corresponds to channel
% label 'c' (due to the selection above). Weirdly, datout now has channels
% a and c, so channel d was removed!
% (issue also holds for trial selection)

end