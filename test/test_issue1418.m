function test_issue1418

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_componentanalysis ft_rejectcomponent

%%

nchan = 10;
ntrl = 20;
ntime = 1000;
fsample = 1000;

data = [];
for i=1:nchan
  data.label{i} = sprintf('eeg%d', i);
end
for i=1:ntrl
  dat = randn(nchan, ntime);
  dat = sign(dat) .* dat.^2; % do something non-linear to make it non-Gaussian
  
  data.time{i} = (1:ntime)/fsample;
  data.trial{i} = dat;
end

nchan = 2;

trig = [];
for i=1:nchan
  trig.label{i} = sprintf('trigger%d', i);
end
for i=1:ntrl
  dat = zeros(nchan, ntime);
  trig.time{i} = (1:ntime)/fsample;
  trig.trial{i} = dat;
end

cfg = [];
data = ft_appenddata(cfg, data, trig);

%%

cfg = [];
cfg.numcomponent = 8;
cfg.method = 'runica';
cfg.channel = {'eeg*'};
% cfg.feedback = 'gui'; % this is to test unrelated improvement that I made with this pull request
comp = ft_componentanalysis(cfg, data);

assert(numel(comp.label)==8);

%%

% this is where some user interaction is actually needed to check that it prints
%   removing 2 components
%   keeping 6 components

cfg = [];
cfg.component = [6 7];
clean1 = ft_rejectcomponent(cfg, comp);
assert(numel(clean1.label)==10);  % this is only the EEG channels

cfg = [];
cfg.component = [6 7];
clean2 = ft_rejectcomponent(cfg, comp, data);
assert(numel(clean2.label)==12);  % this also includes the trigger channels

