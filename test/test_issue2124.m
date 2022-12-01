function test_issue2124

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_rejectartifact

% reported issue is that ft_rejectartifact (with data in the input) fails
% on overlapping data segments, which used to work in the past

% create a minimal data example
data = [];
data.trial{1} = randn(1,100);
data.trial{2} = randn(1,100);
data.trial{3} = randn(1,100);
data.time{1}  = (0:99)./100;
data.time{2}  = (0:99)./100;
data.time{3}  = (0:99)./100;
data.label    = {'label'};
data.sampleinfo = [1 100;51 150;81 180];
data.trialinfo = [1 2 3]';

cfg = [];
cfg.artfctdef.reject   = 'complete';
cfg.artfctdef.zvalue.artifact = [11 20];
dataout = ft_rejectartifact(cfg, data);

cfg = [];
cfg.artfctdef.reject   = 'partial';
cfg.artfctdef.zvalue.artifact = [11 20];
cfg.artfctdef.minaccepttim = 0.05;
dataout = ft_rejectartifact(cfg, data);
assert(numel(dataout.trial)==4);

cfg = [];
cfg.artfctdef.reject   = 'partial';
cfg.artfctdef.zvalue.artifact = [11 20;111 120];
cfg.artfctdef.minaccepttim = 0.05;
dataout = ft_rejectartifact(cfg, data);
assert(numel(dataout.trial)==6);

