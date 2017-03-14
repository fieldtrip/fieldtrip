function test_ft_topoplotCC

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_topoplotCC

% this function tests the ft_topoplotCC function and should display a
% figure with the ctf151 layout and a few arrows.  It does no explicit
% comparisons, but only checks whether it runs through

% create artificial data
conn.label = {
  'MLT22'
  'MRT22'
  'MZO01'
  };
conn.dimord = 'chan_chan_freq';
conn.freq = 1:10;
conn.cohspctrm = rand(3,3, 10);

% make a plot
cfg = [];
cfg.layout = 'CTF151.lay';
cfg.layout = ft_prepare_layout(cfg);

cfg.foi = 1;
cfg.alphaparam = 'cohspctrm';

cfg.arrowhead   = 'stop';
cfg.arrowoffset = 0.03;
cfg.arrowlength = 0.8;

cfg.newfigure   = 'yes';
cfg.alphaparam  = 'cohspctrm';
ft_topoplotCC(cfg, conn);

cfg.arrowhead   = 'none';
cfg.arrowoffset = 0;
cfg.arrowlength = 0.8;

cfg.newfigure   = 'no';
cfg.alphaparam  = 'cohspctrm';
ft_topoplotCC(cfg, conn);
