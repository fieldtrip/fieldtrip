function test_ft_sourceinterpolate

% MEM 4500mb
% WALLTIME 00:10:00

% TEST ft_sourceinterpolate ft_sourceplot

% See also test_bug2769 which goes over more interpolation options using fake data

clear all
close all

%%

srf_lo = ft_read_headshape('cortex_5124.surf.gii');
srf_hi = ft_read_headshape('cortex_20484.surf.gii');
tmp = load('standard_sourcemodel3d10mm.mat'); vol_lo = tmp.sourcemodel;
tmp = load('standard_sourcemodel3d4mm.mat');  vol_hi = tmp.sourcemodel;

figure
ft_plot_mesh(srf_lo);
figure
ft_plot_mesh(srf_hi);
figure
ft_plot_mesh(vol_lo.pos);
figure
ft_plot_mesh(vol_hi.pos);

%% the functional value at each dipole location is a scalar, i.e. 0 dimensions

source0d_vol_lo = [];
source0d_vol_lo.pos = vol_lo.pos;
source0d_vol_lo.dim = vol_lo.dim;
source0d_vol_lo.pow = vol_lo.pos(:,3);
source0d_vol_lo.powdimord = 'pos';

source0d_vol_hi = [];
source0d_vol_hi.pos = vol_hi.pos;
source0d_vol_hi.dim = vol_hi.dim;
source0d_vol_hi.pow = vol_hi.pos(:,3);
source0d_vol_hi.powdimord = 'pos';

source0d_srf_lo = [];
source0d_srf_lo.pos = srf_lo.pos;
source0d_srf_lo.tri = srf_lo.tri;
source0d_srf_lo.pow = srf_lo.pos(:,3);
source0d_srf_lo.powdimord = 'pos';

source0d_srf_hi = [];
source0d_srf_hi.pos = srf_hi.pos;
source0d_srf_hi.tri = srf_hi.tri;
source0d_srf_hi.pow = srf_hi.pos(:,3);
source0d_srf_hi.powdimord = 'pos';

%%

cfg = [];
cfg.parameter = 'pow';
interp0d_vol2vol = ft_sourceinterpolate(cfg, source0d_vol_lo, source0d_vol_hi);
interp0d_vol2srf = ft_sourceinterpolate(cfg, source0d_vol_lo, source0d_srf_hi);
interp0d_srf2vol = ft_sourceinterpolate(cfg, source0d_srf_lo, source0d_vol_hi);
interp0d_srf2srf = ft_sourceinterpolate(cfg, source0d_srf_lo, source0d_srf_hi);

%%

close all
cfg = [];
cfg.funparameter = 'pow';
% these only work on volume data
cfg.method = 'ortho';
ft_sourceplot(cfg, interp0d_vol2vol);
cfg.method = 'slice';
ft_sourceplot(cfg, interp0d_vol2vol);
cfg.method = 'glassbrain'; 
ft_sourceplot(cfg, interp0d_vol2vol);
cfg.method = 'ortho';
ft_sourceplot(cfg, interp0d_srf2vol);
cfg.method = 'slice';
ft_sourceplot(cfg, interp0d_srf2vol);
cfg.method = 'glassbrain'; 
ft_sourceplot(cfg, interp0d_srf2vol);

% these only work on surface data
cfg.method = 'surface';
ft_sourceplot(cfg, interp0d_vol2srf);
cfg.method = 'surface';
ft_sourceplot(cfg, interp0d_srf2srf);

% this works on both, but is rather slow
cfg.method = 'vertex';
ft_sourceplot(cfg, interp0d_srf2vol);
ft_sourceplot(cfg, interp0d_vol2srf);


%% the functional value at each dipole location is a vector, i.e. 1 dimensions

source1d_vol_lo = [];
source1d_vol_lo.pos = vol_lo.pos;
source1d_vol_lo.dim = vol_lo.dim;
source1d_vol_lo.time = 1:7;
for t=1:7
  source1d_vol_lo.pow(:,t) = t*vol_lo.pos(:,3);
end
source1d_vol_lo.powdimord = 'pos_time';


source1d_vol_hi = [];
source1d_vol_hi.pos = vol_hi.pos;
source1d_vol_hi.dim = vol_hi.dim;
source1d_vol_hi.time = 1:7;
for t=1:7
  source1d_vol_hi.pow(:,t) = t*vol_hi.pos(:,3);
end
source1d_vol_hi.powdimord = 'pos_time';

source1d_srf_lo = [];
source1d_srf_lo.pos = srf_lo.pos;
source1d_srf_lo.tri = srf_lo.tri;
source1d_srf_lo.time = 1:7;
for t=1:7
  source1d_srf_lo.pow(:,t) = t*srf_lo.pos(:,3);
end
source1d_srf_lo.powdimord = 'pos_time';

source1d_srf_hi = [];
source1d_srf_hi.pos = srf_hi.pos;
source1d_srf_hi.tri = srf_hi.tri;
source1d_srf_hi.time = 1:7;
for t=1:7
  source1d_srf_hi.pow(:,t) = t*srf_hi.pos(:,3);
end
source1d_srf_hi.powdimord = 'pos_time';

%%

cfg = [];
cfg.parameter = 'pow';
interp1d_vol2vol = ft_sourceinterpolate(cfg, source1d_vol_lo, source1d_vol_hi);
interp1d_vol2srf = ft_sourceinterpolate(cfg, source1d_vol_lo, source1d_srf_hi);
interp1d_srf2vol = ft_sourceinterpolate(cfg, source1d_srf_lo, source1d_vol_hi);
interp1d_srf2srf = ft_sourceinterpolate(cfg, source1d_srf_lo, source1d_srf_hi);

%%

close all
cfg = [];
cfg.funparameter = 'pow';
cfg.method = 'ortho';
ft_sourceplot(cfg, interp1d_vol2vol);
%cfg.method = 'surface'; %-> this does not work, but is not due to ft_sourceinterpolate
%ft_sourceplot(cfg, interp1d_vol2srf);
cfg.method = 'ortho';
ft_sourceplot(cfg, interp1d_srf2vol);
%cfg.method = 'surface';
%ft_sourceplot(cfg, interp1d_srf2srf);


%% the functional value at each dipole location is a matrix, i.e. 2 dimensions

source2d_vol_lo = [];
source2d_vol_lo.pos = vol_lo.pos;
source2d_vol_lo.dim = vol_lo.dim;
source2d_vol_lo.freq = 1:5;
source2d_vol_lo.time = 1:7;
for f=1:5
  for t=1:7
    source2d_vol_lo.pow(:,f,t) = f*t*vol_lo.pos(:,3);
  end
end
source2d_vol_lo.powdimord = 'pos_freq_time';


source2d_vol_hi = [];
source2d_vol_hi.pos = vol_hi.pos;
source2d_vol_hi.dim = vol_hi.dim;
source2d_vol_hi.freq = 1:5;
source2d_vol_hi.time = 1:7;
for f=1:5
  for t=1:7
    source2d_vol_hi.pow(:,f,t) = f*t*vol_hi.pos(:,3);
  end
end
source2d_vol_hi.powdimord = 'pos_freq_time';

source2d_srf_lo = [];
source2d_srf_lo.pos = srf_lo.pos;
source2d_srf_lo.tri = srf_lo.tri;
source2d_srf_lo.freq = 1:5;
source2d_srf_lo.time = 1:7;
for f=1:5
  for t=1:7
    source2d_srf_lo.pow(:,f,t) = f*t*srf_lo.pos(:,3);
  end
end
source2d_srf_lo.powdimord = 'pos_freq_time';

source2d_srf_hi = [];
source2d_srf_hi.pos = srf_hi.pos;
source2d_srf_hi.tri = srf_hi.tri;
source2d_srf_hi.freq = 1:5;
source2d_srf_hi.time = 1:7;
for f=1:5
  for t=1:7
    source2d_srf_hi.pow(:,f,t) = f*t*srf_hi.pos(:,3);
  end
end
source2d_srf_hi.powdimord = 'pos_freq_time';

%%

cfg = [];
cfg.parameter = 'pow';
interp2d_vol2vol = ft_sourceinterpolate(cfg, source2d_vol_lo, source2d_vol_hi);
interp2d_vol2srf = ft_sourceinterpolate(cfg, source2d_vol_lo, source2d_srf_hi);
interp2d_sfr2vol = ft_sourceinterpolate(cfg, source2d_srf_lo, source2d_vol_hi);
interp2d_srf2srf = ft_sourceinterpolate(cfg, source2d_srf_lo, source2d_srf_hi);

%%

close all
cfg = [];
cfg.funparameter = 'pow';
cfg.method = 'ortho';
ft_sourceplot(cfg, interp2d_vol2vol);
%cfg.method = 'surface';
%ft_sourceplot(cfg, interp2d_vol2srf);
cfg.method = 'ortho';
ft_sourceplot(cfg, interp2d_sfr2vol);
%cfg.method = 'surface';
%ft_sourceplot(cfg, interp2d_srf2srf);

%%

[sphere_lo.pnt, sphere_lo.tri] = icosahedron162;
[sphere_hi.pnt, sphere_hi.tri] = icosahedron642;

source0d_sphere_lo = [];
source0d_sphere_lo.pos = sphere_lo.pnt;
source0d_sphere_lo.tri = sphere_lo.tri;
source0d_sphere_lo.pow = sphere_lo.pnt(:,3);
source0d_sphere_lo.powdimord = 'pos';

source1d_sphere_lo = [];
source1d_sphere_lo.pos = sphere_lo.pnt;
source1d_sphere_lo.tri = sphere_lo.tri;
source1d_sphere_lo.time = 1:7;
for t=1:7
  source1d_sphere_lo.pow(:,t) = t*sphere_lo.pnt(:,3);
end
source1d_sphere_lo.powdimord = 'pos_time';

source2d_sphere_lo = [];
source2d_sphere_lo.pos = sphere_lo.pnt;
source2d_sphere_lo.tri = sphere_lo.tri;
source2d_sphere_lo.freq = 1:5;
source2d_sphere_lo.time = 1:7;
for f=1:5
  for t=1:7
    source2d_sphere_lo.pow(:,f,t) = f*t*sphere_lo.pnt(:,3);
  end
end
source2d_sphere_lo.powdimord = 'pos_freq_time';

source0d_sphere_hi = [];
source0d_sphere_hi.pos = sphere_lo.pnt;
source0d_sphere_hi.tri = sphere_lo.tri;
source0d_sphere_hi.orig.pos = sphere_hi.pnt;
source0d_sphere_hi.orig.tri = sphere_hi.tri;
% source0d_sphere_hi.pow = sphere_hi.pnt(:,3);
% source0d_sphere_hi.powdimord = 'pos';

cfg = [];
cfg.parameter = 'pow';
cfg.method = 'smudge';
interp0d_smudge = ft_sourceinterpolate(cfg, source0d_sphere_lo, source0d_sphere_hi);
interp1d_smudge = ft_sourceinterpolate(cfg, source1d_sphere_lo, source0d_sphere_hi);
interp2d_smudge = ft_sourceinterpolate(cfg, source2d_sphere_lo, source0d_sphere_hi);

close all
cfg = [];
cfg.funparameter = 'pow';
cfg.method = 'surface';
ft_sourceplot(cfg, interp0d_smudge);
%ft_sourceplot(cfg, interp1d_smudge);
%ft_sourceplot(cfg, interp2d_smudge);

