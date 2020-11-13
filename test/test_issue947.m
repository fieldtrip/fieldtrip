function test_issue947

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

% the warning state is kept in the global ft_default structure
global ft_default


%%

ft_default = [];

% both should show
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

%%

ft_warning('off', 'a:a')

% only b should show
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

%%

ft_default = [];

ft_warning('off')

% none of them should show
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

%%

ft_warning('on')

% both of them should show again
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

%%

ft_default = [];

ft_warning('off', 'a:a')
ft_warning('off', 'b:b')

% none of them should show
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

%%

ft_default = [];

wsa = ft_warning('off', 'a:a')
wsb = ft_warning('off', 'b:b')

% none of them should show
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

% revert the state of a
ft_warning(wsa)

% only a should show
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

% revert the state of b
ft_warning(wsb)

% both should show
ft_warning('a:a', 'a')
ft_warning('b:b', 'b')

%%

ft_default = [];

% the desired behavior in the following is not really well defined

wsa = ft_warning('off', 'a:a')
ft_warning('off')
ft_warning(wsa) % revert the state for a
ft_warning('a:a', 'a') % this shows
ft_warning('b:b', 'b') % this not
ft_warning('c:c', 'c') % this not

