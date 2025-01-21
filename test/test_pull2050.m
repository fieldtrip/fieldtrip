function test_pull2050

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_event
% DATA no

filename = dccnpath(fullfile('/project/3031000.02/test/issue2047', 'iCM_003_20181205_081017_Cz-epo.fif'));
event    = ft_read_event(filename); % this should not fail
