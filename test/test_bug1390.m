function test_bug1390

% TEST test_bug1390
% TEST ft_timelockanalysis ft_datatype_raw

load /home/common/matlab/fieldtrip/data/test/bug1390.mat

% The following caused a problem in versions prior to 27 March 2012
% that expressed itself as
%   ??? Error using ==> ft_timelockanalysis at 199
%   data has variable trial lengths, you specified not to accept that !
timelock = ft_timelockanalysis(cfg, data);
