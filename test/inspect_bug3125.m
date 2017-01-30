function insspect_bug3125

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3125'))

load('bnd4_1922_corrected_fiducials.mat')
load('bnd4_1922_corrected.mat')

elec1005 = ft_read_sens('standard_1005.elc');

%%

headshape.pos = [];
headshape.fid.pos = [
  NAS
  LPA
  RPA
  ];
headshape.fid.label = {
  'NAS'
  'LPA'
  'RPA'
  };


figure;
ft_plot_sens(elec1005)

figure;
ft_plot_headshape(headshape)
ft_plot_mesh(new_bnd);


%%

cfg = [];
cfg.method = 'fiducial';
cfg.target.elecpos(1,:) = NAS;
cfg.target.elecpos(2,:) = LPA;
cfg.target.elecpos(3,:) = RPA;
cfg.target.chanpos = cfg.target.elecpos;
cfg.target.label = {'Nz','LPA','RPA'};
cfg.target.unit = 'm';
cfg.elec = elec1005;

elec = ft_electroderealign(cfg);
elec = ft_convert_units(elec,'m');

figure
ft_plot_mesh(new_bnd(1))
ft_plot_sens(elec)

%%

cfg = [];
cfg.elec = elec;
cfg.method = 'project';
cfg.headshape = new_bnd(1);
elec = ft_electroderealign(cfg, elec);

figure
ft_plot_mesh(new_bnd(1))
ft_plot_sens(elec)

