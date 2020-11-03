function test_issue1568

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_databrowser

global ft_default
ft_default.representation = 'table';


%%

ntrials = 10;
nchans = 16;
fsample = 1000;
trialgap = 0;

data = [];
for i=1:nchans
  data.label{i} = sprintf('chan%d', i);
end
for i=1:ntrials
  data.time{i} = (1:fsample)/fsample;
  data.trial{i} = randn(nchans, fsample);
  data.sampleinfo(i,1) = (i-1)*fsample + 1 + (i*trialgap); % create some gaps between the trials
  data.sampleinfo(i,2) = data.sampleinfo(i,1) + fsample - 1;
end

if strcmp(ft_default.representation, 'numeric')
  data.trialinfo = [1:ntrials]';
else
  trialnum = arrayfun(@num2str, 1:ntrials, 'UniformOutput', false)';
  data.trialinfo = table(trialnum);
end

%%
close all

cfg = [];
cfg.ylim = [-6 6];
cfg.viewmode = 'vertical';

tmp1.begsample = 101;
tmp1.endsample = 200;
tmp1.channel = 'chan1';
tmp1.value = 1;

tmp2.begsample = 301;
tmp2.endsample = 400;
tmp2.channel = 'chan2';
tmp2.value = 2;

tmp3.begsample = 501;
tmp3.endsample = 600;
tmp3.channel = 'chan3';
tmp3.value = 3;

% let's use three different types of artifacts, just like visual, eog, muscle, clip, etc.
cfg.artfctdef.art1.artifact = struct2table(tmp1);
cfg.artfctdef.art2.artifact = struct2table(tmp2);
cfg.artfctdef.art3.artifact = struct2table(tmp3);

event(1).type     = 'ttl1';
event(1).value    = 1;
event(1).sample   = 701;
event(1).duration = 10;
event(1).offset   = 0;

event(2).type     = 'ttl2';
event(2).value    = 2;
event(2).sample   = 721;
event(2).duration = 10;
event(2).offset   = 0;

event(3).type     = 'ttl3';
event(3).value    = 3;
event(3).sample   = 741;
event(3).duration = 10;
event(3).offset   = 0;

cfg.event = event;

% cfg.ploteventlabels = 'type=value';
% ft_databrowser(cfg, data);
% 
% cfg.ploteventlabels = 'type';
% ft_databrowser(cfg, data);

cfg.ploteventlabels = 'value';
ft_databrowser(cfg, data);

