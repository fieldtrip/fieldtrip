function inspect_issue1674

% WALLTIME 00:10:00
% MEM 2gb

%%

elec = [];
elec.label = {
  'Fpz'
  'T7'
  'Cz'
  'T8'
  'Oz'
  };
elec.elecpos = [
  80   0   0 % a
   0  80   0 % l
   0   0  80 % s
   0 -80   0
 -80   0   0
  ];
elec.unit = 'mm';
elec.coordsys = 'eeglab';

%%

cfg = [];
cfg.elec = elec;
layout = ft_prepare_layout(cfg);
ft_plot_layout(layout);

