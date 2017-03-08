function test_bug1390

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_timelockanalysis ft_datatype_raw

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1390.mat'));

% The following caused a problem in versions prior to 27 March 2012
% that expressed itself as
%   ??? Error using ==> ft_timelockanalysis at 199
%   data has variable trial lengths, you specified not to accept that !
timelock = ft_timelockanalysis(cfg, data);

% this is the follow up in comment 8
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1390c8.mat'));
timelock = ft_timelockanalysis(cfg, data);

