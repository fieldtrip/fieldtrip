function test_bug367

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_megrealign ft_read_sens ft_dipolesimulation ft_timelockanalysis ft_multiplotER

% this test script is for http://bugzilla.fcdonders.nl/show_bug.cgi?id=367
% and was constructed from http://www.fieldtriptoolbox.org/example/megrealign

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% the two example files with the 151 and 275 channel CTF specifications are available from
% ftp://ftp.fcdonders.nl/pub/fieldtrip/example/megrealign

grad151 = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/bug367/ctf151.mat'));
grad275 = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/bug367/ctf275.mat'));

vol = [];
vol.r = 12;
vol.o = [0 0 4];
vol.unit = 'cm';

% in this example the dipole is carefully positioned on the same depth as
% where the dipole layer for the interpolation will be located. The center
% of the head is at [0 0 4], the radius is 12 cm, and the dipole layer that is
% specified in ft_megrealign is 2.5 cm shifted inward from the head surface
cfg = [];
cfg.dip.pos = [0 0 13.5];  % 4 + 12 - 2.5
cfg.dip.frequency = 1;
cfg.vol = vol;
cfg.grad = grad151;
data151 = ft_dipolesimulation(cfg);
cfg.grad = grad275;
data275 = ft_dipolesimulation(cfg);

cfg = [];
avg151 = ft_timelockanalysis(cfg, data151);
avg275 = ft_timelockanalysis(cfg, data275);

% here I'll not only realign the 151 channel MEG data to the 275 channel and
% vice versa, but also apply the ft_megrealign function from 151 to 151. In the
% realignment the signal will be slightly distorted as a comparison of avg151_151
% to avg151 would show.
cfg = [];
cfg.inwardshift = 3;
cfg.vol = vol;
cfg.template{1} = grad151; avg151_151 = ft_timelockanalysis([], ft_megrealign(cfg, avg151));
cfg.template{1} = grad275; avg151_275 = ft_timelockanalysis([], ft_megrealign(cfg, avg151));
cfg.template{1} = grad151; avg275_151 = ft_timelockanalysis([], ft_megrealign(cfg, avg275));
cfg.template{1} = grad275; avg275_275 = ft_timelockanalysis([], ft_megrealign(cfg, avg275));

% plot the realigned datasets
cfg = [];
figure; ft_multiplotER(cfg, avg151_151);
title('Interpolated from 151 to 151');

figure; ft_multiplotER(cfg, avg151_275);
title('Interpolated from 151 to 151');

figure; ft_multiplotER(cfg, avg275_151);
title('Interpolated from 151 to 151');

figure; ft_multiplotER(cfg, avg275_275);
title('Interpolated from 151 to 151');

% now plot them together
cfg = [];
figure; ft_multiplotER(cfg, avg151, avg151_151, avg275_151);
figure; ft_multiplotER(cfg, avg275, avg275_275, avg151_275);


