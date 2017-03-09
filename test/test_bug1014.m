function test_bug1014

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_checkdata ft_appendtimelock ft_selectdata

a.time   = [1 2];
a.label  = {'chan'};
a.dimord = 'chan_time';
a.avg    = [1 1];
a.trial  = reshape([3 3; 0 0; 0 0], [3 1 2]);

b.time   = [1 2];
b.label  = {'chan'};
b.dimord = 'chan_time';
b.avg    = [2 2];
b.trial  = reshape([6 6; 0 0; 0 0], [3 1 2]);

% the following is to work around
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1013

a.sampleinfo = [1 2; 3 4; 5 6];
b.sampleinfo = [1 2; 3 4; 5 6];

% re-test bug 1013
ft_checkdata(a, 'datatype', 'timelock', 'hassampleinfo', 'ifmakessense');
ft_checkdata(b, 'datatype', 'timelock', 'hassampleinfo', 'ifmakessense');

c = ft_appendtimelock([], a, b);

if isfield(c, 'avg') && ~isfield(c, 'trial')
  error('the result should have trial rather than avg');
end

a = rmfield(a, 'trial');
b = rmfield(b, 'trial');
c = ft_appendtimelock([], a, b);

if isfield(c, 'avg')
  error('the result should not have an average');
end
