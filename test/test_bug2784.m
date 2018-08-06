function test_bug2784

% MEM 150mb
% WALLTIME 00:10:00

% TEST ft_mvaranalysis

data = [];
data.trial{1} = randn(2,100);
data.trial{2} = randn(2,101);
data.time{1}  = (0:99)./100+0.001;
data.time{2}  = (0:100)./100;
data.label    = {'chan01';'chan02'};

try
  cfg = [];
  ft_mvaranalysis(cfg,data);
catch me
  if strcmp(me.message, 'time axes of all trials should be identical')
    fprintf('error caught by function\n');
  end
end

cfg.t_ftimwin = 0.5;
cfg.toi       = 0.25;
ft_mvaranalysis(cfg,data);

% load in the data provided by Tyler Grummett
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2784.mat');
load(filename);

cfg = [];
mdata = ft_mvaranalysis(cfg, temp);

% the issue Tyler had could be reproduced and was caused by a combination
% of things: the tfwin being 501 samples long, as well as the time axes of
% the trials being different. ft_mvaranalysis has been adjusted to more
% robustly deal with these cases.
