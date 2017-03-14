function test_bug2377

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_datatype_sens ft_compute_leadfield

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2377/eeg_lf_scaling.mat'));

sens = rmfield(sens, 'tra');

% initially the electrodes are not on the skin surface
figure
ft_plot_vol(vol);
ft_plot_sens(sens);
camlight
alpha 0.5

[vol, sens] = ft_prepare_vol_sens(vol, sens);

% now they are on the skin surface
figure
ft_plot_vol(vol);
ft_plot_sens(sens);
camlight
alpha 0.5

assert(strcmp(vol.unit, 'm'));
assert(strcmp(sens.unit, 'm'));
assert(all(strcmp(sens.chanunit, 'V')));

%% compute some leadfields

lf1 = ft_compute_leadfield([0 0 0.05], sens, vol, 'chanunit', [], 'dipoleunit', []); % default units, i.e. SI
lf2 = ft_compute_leadfield([0 0 0.05], sens, vol, 'chanunit', repmat({'uV'}, 1, 128), 'dipoleunit', []);
lf3 = ft_compute_leadfield([0 0 0.05], sens, vol, 'chanunit', [], 'dipoleunit', 'nA*m');
lf4 = ft_compute_leadfield([0 0 0.05], sens, vol, 'chanunit', repmat({'uV'}, 1, 128), 'dipoleunit', 'nA*m');

figure
% ft_plot_vol(vol);
ft_plot_topo3d(sens.chanpos, lf4(:,1));
colorbar


%%  check the relative amplitudes

n1 = norm(lf1);
n2 = norm(lf2) % 1e6  times n1;
n3 = norm(lf3) % 1e-9 times n1;
n4 = norm(lf4) % 1e-3 times n1;

assert(abs(n2/n1-1e6)/1e6   < 100*eps); 
assert(abs(n3/n1-1e-9)/1e-9 < 100*eps); 
assert(abs(n4/n1-1e-3)/1e-3 < 100*eps); 

%% FIXME this does not mean that the absolte values are correct
