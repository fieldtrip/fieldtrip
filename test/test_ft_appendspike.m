function test_ft_appendspike

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_appendspike

spike = [];
spike.label = {'unit1'  'unit2'  'unit3'};
spike.timestamp = {linspace(0,100,100) linspace(0,200,200) linspace(0,300,300)};
spike.trial = {randn(1,100) randn(1,200) randn(1,300)};
spike.trialtime = [linspace(0,30,300)', linspace(0,30,300)'+0.5];

spike2 = spike;
spike.label = {'unit12'  'unit22'  'unit32'};
spike2.trial = {randn(1,100) randn(1,200) randn(1,300)};

cfg = [];
cfg.trl = spike.trialtime;
spikeout = ft_appendspike(cfg,spike,spike2);