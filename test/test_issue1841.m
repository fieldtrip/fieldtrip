function test_issue1841

% WALLTIME 00:30:00
% MEM 4gb
% DEPENDENCY ft_scalpcurrentdensity ft_channelrepair

%%
[ftver, ftdir] = ft_version;
elec = ft_read_sens(fullfile(ftdir, 'template/electrode/standard_1020.elc'));

cd(fullfile(ftdir, 'private')); % this needs to be done, otherwise private functions cannot be tested

order  = 4;
lambda = 1e-5;
V1     = randn(size(elec.chanpos,1),1);

[V2, L2, L1]           = splint(elec.chanpos, V1, elec.chanpos, order, 500, lambda);
[W, Gss, gx2, hx2]     = sphericalSplineInterpolate(elec.chanpos', elec.chanpos', lambda, order, 'slap');

figure;plot(L1, W*V1, 'o'); % this should be the same, or not?

