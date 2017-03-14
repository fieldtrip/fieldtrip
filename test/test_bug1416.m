function test_bug1416

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_setopt

opt = {'opt1' 1 'opt2' 2};

optnew = ft_setopt(opt, 'opt2', 4);
assert(isequal(optnew{1}, 'opt1'));
assert(isequal(optnew{2}, 1));
assert(isequal(optnew{3}, 'opt2'));
assert(isequal(optnew{4}, 4));

optnew = ft_setopt(opt, 'opt3', 3);
assert(isequal(optnew{1}, 'opt1'));
assert(isequal(optnew{2}, 1));
assert(isequal(optnew{3}, 'opt2'));
assert(isequal(optnew{4}, 2));
assert(isequal(optnew{5}, 'opt3'));
assert(isequal(optnew{6}, 3));

opt = {'opt1' 1 'opt1' 2}; % this is an ambigous (invalid?) specification

optnew = ft_setopt(opt, 'opt1', 3); % this should remove the two old ones and add a single new one
assert(isequal(optnew{1}, 'opt1'));
assert(isequal(optnew{2}, 3));
assert(length(optnew)==2);
