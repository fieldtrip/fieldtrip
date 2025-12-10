function test_issue1425

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preprocessing ft_read_data
% DATA private

filename = dccnpath('/project/3031000.02/test/issue1425/X_X_e515c5ac-6301-4acd-8a69-fb208d5fd097_0014.edf');

hdr   = ft_read_header(filename);
event = ft_read_event(filename);

% the third event is an 'XLSpike', which, in line with the EDF browser,
% should have a timestamp of about 1.951, and not 2.345
% NOTE that as per the reported issue (1425) the EDF browser returns
% 1.9511, whereas FT seems to return 1.9511719
assert(isequal(round(event(3).timestamp*10000), 19512) && isequal(event(3).value, 'XLSpike'));