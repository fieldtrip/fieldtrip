function test_issue1841

% WALLTIME 00:30:00
% MEM 4gb
% DEPENDENCY ft_scalpcurrentdensity ft_channelrepair

%%
% Loads the 10-20 sensor definition.
[ftver, ftdir] = ft_version;
elec = ft_read_sens(fullfile(ftdir, 'template/electrode/standard_1020.elc'));

cd(fullfile(ftdir, 'private')); % this needs to be done, otherwise private functions cannot be tested

% Centers the sensor definition.
center = fitsphere ( elec.chanpos );
elec.chanpos = elec.chanpos - center;

% Sets the parameters.
order  = 4;
lambda = 1e-5;
Vi     = randn(size(elec.chanpos,1),1);

% Runs all the functions.
[Vo, Lo, ~]   = splint(elec.chanpos, Vi, elec.chanpos, order, 500, lambda);
[WV, ~, ~, ~] = sphericalSplineInterpolate(elec.chanpos', elec.chanpos', lambda, order, 'spline');
[WL, ~, ~, ~] = sphericalSplineInterpolate(elec.chanpos', elec.chanpos', lambda, order, 'slap');
[WVo, WLo]    = sphsplint(elec.chanpos, elec.chanpos, order, 500, lambda);

% Compares the output of the original functions.
figure
subplot 121
plot(Vo, WV*Vi, 'o'),title('Potential')
subplot 122
plot(Lo, WL*Vi, 'o'),title('Laplacian')

% Compares with the output of the new function.
figure
subplot 121
plot(Vo, WVo*Vi, 'o'),title('Potential')
subplot 122
plot(Lo, WLo*Vi, 'o'),title('Laplacian')
