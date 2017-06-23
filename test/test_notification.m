function test_notification

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_debug ft_info ft_notice ft_warning ft_notification

subfunction(@ft_debug)
subfunction(@ft_info)
subfunction(@ft_notice)
subfunction(@ft_warning)

% ft_error cannot be tested the same way, since the error stops everything
errortest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function errortest

try
  ft_error('off');        % this should not make a difference
  ft_error('something');  % note that the line number is used further down
  ok = false;
catch
  me = ft_error('last')
  assert(strcmp(me.message, 'something'))
  ok = true;
end

if ~ok
  error('ft_error failed');
else
  fprintf('========================================\n')
  fprintf('ft_error looks good\n');
  fprintf('========================================\n')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subfunction(notice)

fprintf('========================================\n')
fprintf('the following should not print anything\n')
fprintf('========================================\n')
notice('off', 'FieldTrip:something');
notice('FieldTrip:something', 'notice %d %d %d\n', 1, 2, 3);

fprintf('========================================\n')
fprintf('the following should print the query status\n')
fprintf('========================================\n')
notice('query', 'FieldTrip:something');

fprintf('========================================\n')
fprintf('the following should print "notice 4, 5, 6"\n')
fprintf('========================================\n')
notice('on', 'FieldTrip:something');
notice('FieldTrip:something', 'notice %d, %d, %d\n', 4, 5, 6);

fprintf('========================================\n')
fprintf('the following should print "notice backtrace" with a backtrace\n')
fprintf('========================================\n')
notice('on', 'backtrace');
notice('FieldTrip:something', 'notice backtrace\n');
notice('off', 'backtrace');

fprintf('========================================\n')
fprintf('the following should print "pause 1", "4"\n')
fprintf('========================================\n')
notice('once', 'FieldTrip:once');
notice('timeout', 3);
for i=1:5
  pause(1);
  notice('FieldTrip:once', 'pause %d\n', i);
end

notice('clear')
notice('on', '1')
notice('on', '2')
notice('off', 'a')
notice('off', 'b')
fprintf('========================================\n')
fprintf('the following should show 1 and 2 as ''on''\n')
fprintf('========================================\n')
notice('off')
notice('query')
fprintf('========================================\n')
fprintf('the following should show a and b as ''off''\n')
fprintf('========================================\n')
notice('on')
notice('query')

notice('clear')
fprintf('========================================\n')
fprintf('the following should show ''a''\n')
fprintf('========================================\n')
notice('a');
notice('last')

