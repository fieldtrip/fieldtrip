function test_old_fixcsd

% MEM 1gb
% WALLTIME 00:10:00


% this script tests the functionality of fixcsd/fixcoh in checkdata

% generate some data
rand('twister', 20070425);
freq        = [];
freq.dimord = 'rpttap_chan_freq';
freq.fourierspctrm = randn(30,2,10) + i.*randn(30,2,10);
freq.freq   = 1:10;
freq.label  = {'a';'b'};
freq.cumtapcnt = 5*ones(6,1);
freq.cfg    = [];

channelcmb = ft_channelcombination({'all' 'all'}, freq.label, 1);

f1 = ft_checkdata(freq, 'cmbrepresentation', 'full');
f2 = ft_checkdata(freq, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', channelcmb(1,:));
f3 = ft_checkdata(freq, 'cmbrepresentation', 'sparse', 'channelcmb', channelcmb);

f1f = f1;
f2f = ft_checkdata(f2, 'cmbrepresentation', 'full');
f3f = ft_checkdata(f3, 'cmbrepresentation', 'full');

f1f2  = ft_checkdata(f1f, 'cmbrepresentation', 'sparse');
f1f2f = ft_checkdata(f1f2, 'cmbrepresentation', 'full');

c3 = f3;
c3.cohspctrm = f3.crsspctrm;
c3f = ft_checkdata(c3, 'cmbrepresentation', 'full');
c3s = ft_checkdata(c3f, 'cmbrepresentation', 'sparse');

c1 = f1f2;
c1.cohspctrm = c1.crsspctrm;
c1f = ft_checkdata(c1, 'cmbrepresentation', 'full');
c1s = ft_checkdata(c1f, 'cmbrepresentation', 'sparse');

cfg = [];
cfg.method = 'coh';
cfg.channelcmb = channelcmb(1,:);
coh = ft_connectivityanalysis(cfg, freq);
coh2 = checkdata(coh, 'cmbrepresentation', 'sparse');

cfg = [];
cfg.method = 'granger';

g = ft_connectivityanalysis(cfg, freq);
g2 = ft_checkdata(g, 'cmbrepresentation', 'sparse');

%this part tests the code which aims at speeding up the fourier->full conversion if all(freq.cumtapcnt==freq.cumtapcnt(1))


% generate some data
rand('twister', 20070425);
freq        = [];
freq.dimord = 'rpttap_chan_freq';
freq.fourierspctrm = randn(9000,30,20) + i.*randn(9000,30,20);
freq.freq   = 1:20;
freq.label  = {'a';'b';'c';'d';'e';'f';'g';'h';'i';'j';'k';'l';'m';'n';'o';'p';'q';'r';'s';'t';'u';'v';'w';'x';'y';'z';'aa';'bb';'cc';'dd'};
freq.cumtapcnt = 3*ones(3000,1);
freq.cfg    = [];

freqx = ft_checkdata(freq, 'cmbrepresentation', 'full');
