function inspect_bug3157

% WALLTIME 0:10:00
% MEM 1gb

%%
ctf151      = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151/Subject01.ds'));
neuromag122 = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag122/jg_single_01raw.fif'));
neuromag306 = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif'), 'senstype', 'meg');

%%

rx = 30;
ry = 30;
rz = 30;

fake = [];
fake.unit = 'mm';
fake.coilpos = [
  10  0 0
  -10 -0 0
  ];
fake.coilpos = ft_warp_apply(rotate([rx ry rz]), fake.coilpos);
fake.coilori = [
  0 0 1
  0 0 1
  ];
fake.coilori = ft_warp_apply(rotate([rx ry rz]), fake.coilori);
fake.label = {
  '1'
  };
fake.tra = [
  1 -1
  ];

figure
ft_plot_sens(fake, 'coil', false, 'unit', 'mm', 'coilsize', 40, 'coilshape', 'square');
ft_plot_sens(fake, 'coil', true,  'unit', 'mm', 'coilsize', 0);
grid on

%%

figure
ft_plot_sens(ctf151, 'coil', false, 'unit', 'mm', 'coilsize', 15, 'coilshape', 'circle', 'chantype', 'meggrad');

%%

figure
ft_plot_sens(neuromag122, 'coil', false, 'unit', 'mm', 'coilsize', 35, 'coilshape', 'square', 'chantype', 'megplanar');

%%

figure
ft_plot_sens(neuromag306, 'coil', false, 'unit', 'mm', 'coilsize', 30, 'coilshape', 'square', 'chantype', 'megplanar');

%%

figure
ft_plot_sens(neuromag306, 'coil', false, 'unit', 'mm', 'coilsize', 30, 'coilshape', 'square', 'chantype', 'megplanar', 'facecolor', 'b', 'facealpha', 0.1);
