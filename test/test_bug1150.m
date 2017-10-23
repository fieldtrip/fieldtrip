function test_bug1150

% MEM 2000mb
% WALLTIME 00:10:00

% TEST ft_sourcestatistics

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1150.mat'));

% the following failed
% stat = ft_sourcestatistics(cfg, temp);
% which could be tracked down to the dimord being incorrect after
% ft_checkdata

% the data seems to be incosistent so the conversion fails.
temp.trial = temp.trial(1:10);
temp.cumtapcnt = temp.cumtapcnt(1:10);
temp.trialinfo = temp.trialinfo(1:10,:);

sourcenew  = ft_checkdata(temp, 'sourcerepresentation', 'new');
if ~isequal(size(sourcenew.pow), [38556 10])
  error('incorrect dimensions');
end

