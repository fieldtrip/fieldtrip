function [badchannels] = detectChannelNaN(cfg,data)

% check out ft_artifact_nan.m
% DETECTCHANNELNAN detects whether one or more channels contain only NaNs
% and if so adds them to cfg.badchannels
badchannels      = ft_getopt(cfg, 'badchannel',     {});

if any(isnan(data.trial{1}(:,1))) % check if any NaNs in first sample of first trial
    fprintf('Found channels starting with NaNs. Will be added to cfg.badchannels\n');
    badChan = ft_channelselection(find(isnan(data.trial{1}(:,1))),data.label);
    for k = 1:length(badChan)
        fprintf('Channel %s contains NaNs, added to bad channels\n',string(badChan(k)));
    end
    badchannels = cat(1, badchannels, badChan);
end
