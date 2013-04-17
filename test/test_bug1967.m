% function test_bug1967

% TEST test_bug1967
% TEST ft_prepare_vol_sens

%% construct a segmentation
nx = 101;
ny = 101;
nz = 101;

seg = zeros(nx, ny, nz);
[posx, posy, posz] = ndgrid(1:nx, 1:ny,1:nz);
posx = posx-(nx-1)/2-1;
posy = posy-(nx-1)/2-1;
posz = posz-(nx-1)/2-1;

distance = sqrt(posx.^2 + posy.^2 + posz.^2);
seg(distance<45) = 1;
seg(distance<40) = 2;
seg(distance<30) = 3;

volume.dim = [nx ny nz];
volume.seg = seg;
volume.seglabel = {'skin' 'skull' 'brain'};
volume.transform = [
  1 0 0 -51 
  0 1 0 -51
  0 0 1 -51 
  0 0 0 1
  ];

cfg = [];
cfg.funparameter = 'seg';
ft_sourceplot(cfg, volume);

%% convert it into a head model
cfg = [];
cfg.tissue = {'skin' 'skull' 'brain'};
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg, volume);

cfg = [];
cfg.tissue = {'skin' 'skull' 'brain'};
cfg.method = 'simbio';
headmodel = ft_prepare_headmodel(cfg, mesh);


