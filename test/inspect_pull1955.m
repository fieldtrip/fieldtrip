function inspect_pull1955

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_databrowser

%%
% create some uniform data

nchan = 10;
ntrial = 30;
nsample = 1000;
fsample = 1000;

data = [];
for i=1:nchan
  data.label{i} = num2str(i);
end

for i=1:ntrial
  data.trial{i} = randn(nchan,nsample);
  data.time{i}  = (1:nsample)/fsample;
end

data.sampleinfo(:,1) = ((1:ntrial)-1)*nsample+1;
data.sampleinfo(:,2) = ((1:ntrial)  )*nsample;

% add an artifact to channel 2, trial 2
data.trial{2}(2,:) = 3 * data.trial{2}(2,:);

%%

cfg = [];
% cfg.channel = 1:5;
% cfg.trials = 1:15;
cfg.method = 'summary';
% cfg.keeptrials = 'no';
cfg.keeptrials = 'no';
% cfg.keeptrials = 'nan';
% cfg.keeptrials = 'zero';
% cfg.keepchannels = 'no';
% cfg.keepchannels = 'yes';
% cfg.keepchannels = 'nan';
% cfg.keepchannels = 'zero';
cfg.keepchannels = 'no';

% neighbours are needed for method=repair
cfg.neighbours(1).label = '1';
cfg.neighbours(1).neighblabel = {'2'};
cfg.neighbours(2).label = '2';
cfg.neighbours(2).neighblabel = {'3'};
cfg.neighbours(3).label = '3';
cfg.neighbours(3).neighblabel = {'4'};
cfg.neighbours(4).label = '4';
cfg.neighbours(4).neighblabel = {'5'};
cfg.neighbours(5).label = '5';
cfg.neighbours(5).neighblabel = {'6'};
cfg.neighbours(6).label = '6';
cfg.neighbours(6).neighblabel = {'7'};
cfg.neighbours(7).label = '7';
cfg.neighbours(7).neighblabel = {'8'};
cfg.neighbours(8).label = '8';
cfg.neighbours(8).neighblabel = {'9'};
cfg.neighbours(9).label = '9';
cfg.neighbours(9).neighblabel = {'10'};
cfg.neighbours(10).label = '10';
cfg.neighbours(10).neighblabel = {'1'};

fprintf('------------------------------------------------\n')

dataclean = ft_rejectvisual(cfg, data)

% these should in all situations match the output data
% these were already present before the pull request
dataclean.cfg.trials
dataclean.cfg.channel

% these should indicate the channels/trials that were marked or initially left out
% only badchannel is new after this pull request
dataclean.cfg.artfctdef.summary.artifact
dataclean.cfg.badchannel



