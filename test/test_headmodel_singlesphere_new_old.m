function test_headmodel_singlesphere_new_old

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_headmodel ft_headmodel_singlesphere ft_prepare_concentricspheres

% generate a unit sphere
[pnt, tri] = icosahedron162;

% create the BEM geometries
geom = [];
geom.bnd(1).pnt = pnt * 100;
geom.bnd(1).tri = tri;
geom.bnd(2).pnt = pnt * 90;
geom.bnd(2).tri = tri;
geom.bnd(3).pnt = pnt * 80;
geom.bnd(3).tri = tri;

elec.chanpos = pnt * 100;
elec.elecpos = pnt * 100;
for i=1:size(pnt,1)
  elec.label{i} = sprintf('%d', i);
end

arg(1).name = 'conductivity';
arg(1).value = {1, 0.33, 0.1};

optarg = constructalloptions(arg);
% random shuffle the configurations
optarg = optarg(randperm(size(optarg,1)), :);

for i=1:size(optarg,1)
  
  arg = optarg(i,:);
  
  % new way - low level: singlesphere
  vol{1} = ft_headmodel_singlesphere(geom.bnd(1),arg{:});
  
  % old way - low level: concentricspheres
  tmpcfg = ft_keyval2cfg(arg{:});
  tmpcfg.headshape = geom.bnd(1);
  vol{2} = ft_prepare_concentricspheres(tmpcfg);
  
  % new way - high level: singlesphere
  tmpcfg = ft_keyval2cfg(arg{:});
  tmpcfg.method = 'singlesphere';
  vol{3} = ft_prepare_headmodel(tmpcfg,geom.bnd(1));
  
  % new way - high level: concentricspheres:
  tmpcfg = ft_keyval2cfg(arg{:});
  tmpcfg.method = 'concentricspheres';
  vol{4} = ft_prepare_headmodel(tmpcfg,geom.bnd(1));
  
  % compute the leadfields for a comparison
  [vol{1}, elec] = ft_prepare_vol_sens(vol{1}, elec);
  [vol{2}, elec] = ft_prepare_vol_sens(vol{2}, elec);
  [vol{3}, elec] = ft_prepare_vol_sens(vol{3}, elec);
  [vol{4}, elec] = ft_prepare_vol_sens(vol{4}, elec);
  lf{1} = ft_compute_leadfield([0 10 60], elec, vol{1});
  lf{2} = ft_compute_leadfield([0 10 60], elec, vol{2});
  lf{3} = ft_compute_leadfield([0 10 60], elec, vol{3});
  lf{4} = ft_compute_leadfield([0 10 60], elec, vol{4});
  
  % compare the leadfields in all possible combinations
  comb = nchoosek(1:numel(vol),2);
  for j=1:size(comb,1)
    chk = comb(j,:);
    err = norm(lf{chk(1)} - lf{chk(2)}) / norm(lf{chk(1)});
    if err>0.001
      error('combination %d %d not successful\n',chk(1),chk(2));
    end
  end
  
end % for different conductivities
