function test_bug3343

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY

%% this used to return ctf
elec = [];
elec.elecpos = randn(3,2);
elec.chanpos = elec.elecpos;
elec.label = {'BR1';'BR2'};

assert(~strcmp(ft_senstype(elec), 'ctf'))
assert( strcmp(ft_senstype(elec), 'eeg')) % because it has elecpos

%% this used to return ctf
elec = [];
elec.elecpos = randn(3,2);
elec.chanpos = elec.elecpos;
elec.label = {'BR1';'BR2'};
elec.chantype = {'eeg';'eeg'};

assert(~strcmp(ft_senstype(elec), 'ctf'))
assert( strcmp(ft_senstype(elec), 'eeg')) % because of the chantype

%% use boolean in detection
elec.elecpos = randn(3,2);
elec.chanpos = elec.elecpos;
elec.label = {'BR1';'BR2'};
elec.chantype = {'seeg';'ecog'};

assert( ft_senstype(elec, 'eeg'))
assert( ft_senstype(elec, 'ieeg'))
assert(~ft_senstype(elec, 'seeg')) % not all channels are like this
assert(~ft_senstype(elec, 'ecog')) % not all channels are like this

%% use boolean in detection
elec.elecpos = randn(3,2);
elec.chanpos = elec.elecpos;
elec.label = {'BR1';'BR2'};
elec.chantype = {'seeg';'seeg'};

assert( ft_senstype(elec, 'eeg'))
assert( ft_senstype(elec, 'ieeg'))
assert( ft_senstype(elec, 'seeg'))
assert(~ft_senstype(elec, 'ecog'))

%% use boolean in detection
elec.elecpos = randn(3,2);
elec.chanpos = elec.elecpos;
elec.label = {'BR1';'BR2'};
elec.chantype = {'ecog';'ecog'};

assert( ft_senstype(elec, 'eeg'))
assert( ft_senstype(elec, 'ieeg'))
assert( ft_senstype(elec, 'ecog'))
assert(~ft_senstype(elec, 'seeg'))
