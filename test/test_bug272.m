function test_bug272

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_timelockanalysis ft_prepare_layout ft_multiplotER ft_topoplotER 

% this script tests bug 272, i.e. incompatibility of ft_multiplotER with
% cfg.inputfile and tests the fix.

% create some ERP-data
cfg = [];
cfg.layout = 'CTF151.lay';
lay = ft_prepare_layout(cfg);

data = [];
data.label    = lay.label(1:151);
data.trial{1} = randn(151,100);
data.time{1}  = (0:99)/100;
data.fsample  = 100;

tname1 = tempname;
cfg = [];
cfg.outputfile = tname1;
tlck = ft_timelockanalysis(cfg, data);

data.trial{1} = randn(151,100);
tname2 = tempname;
cfg.outputfile = tname2;
tlck2 = ft_timelockanalysis(cfg, data);

cfg           = [];
cfg.layout    = lay;
cfg.inputfile = {tname1 tname2};
ft_multiplotER(cfg);

cfg.inputfile = {tname1};
ft_topoplotER(cfg);
cfg.inputfile = {tname1 tname2};
ft_topoplotER(cfg); % should give an error

delete([tname1,'.mat']);
delete([tname2,'.mat']);

