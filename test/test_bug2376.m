function test_bug2376

% MEM 1500mb
% WALLTIME 00:10:00

% TEST channelposition ft_plot_sens ft_plot_vol ft_plot_mesh ft_prepare_vol_sens

[pnt, tri] = icosahedron162;

% create volume conductor models
vol = [];
vol.o = [0 0 0];
vol.r = 10;
vol.unit = 'cm';
vol.type = 'singlesphere';

% create sensor array, note that it is 11 cm instead of 10
sens.pnt = pnt.*11;
sens.pnt(pnt(:,3)<0,:) = [];
for k = 1:size(sens.pnt,1)
  sens.label{k,1} = num2str(k,'%03d');
end
sens.unit = 'cm';

sens = ft_datatype_sens(sens);
assert(isequal(sens.chanpos, sens.elecpos))

% the electrodes should be around the sphere at 1 cm distance
figure
ft_plot_vol(vol); camlight
ft_plot_sens(sens);

[volp, sensp] = ft_prepare_vol_sens(vol, sens);
assert(isequal(sensp.chanpos, sensp.elecpos))

% the electrodes should be right on the sphere
figure
ft_plot_vol(volp); camlight
ft_plot_sens(sensp);
