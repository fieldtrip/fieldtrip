function test_bug670

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_channelselection ft_channelcombination

datachannel = {'a1', 'a2', 'a3', 'b1', 'b2', 'b3'};

% test wether it returns the expected selections
assert(length(ft_channelselection('all', datachannel))==6);
assert(length(ft_channelselection('a1', datachannel))==1);
assert(length(ft_channelselection({'a1', 'a2'}, datachannel))==2);
assert(length(ft_channelselection({'all', 'a2'}, datachannel))==6);
assert(length(ft_channelselection('a*', datachannel))==3);
assert(length(ft_channelselection({'a*', 'b1'}, datachannel))==4);
assert(length(ft_channelselection({'all', '-b1'}, datachannel))==5);
assert(length(ft_channelselection({'all', '-b*'}, datachannel))==3);
assert(length(ft_channelselection({'all', '-*1'}, datachannel))==4);
assert(length(ft_channelselection('*1', datachannel))==2);

% numeric selections should also work
assert(length(ft_channelselection(1, datachannel))==1);
assert(length(ft_channelselection(0, datachannel))==0);
assert(length(ft_channelselection([1 4], datachannel))==2);
assert(length(ft_channelselection(1:6, datachannel))==6);
assert(length(ft_channelselection(1:7, datachannel))==6);

% and also do some tests for channel combinations
assert(size(ft_channelcombination({'a*', 'b*'}, datachannel),1)==9);
assert(size(ft_channelcombination({'all', 'b*'}, datachannel),1)==(5+4+3));
assert(size(ft_channelcombination({'all', 'all'}, datachannel),1)==(6*5/2));

