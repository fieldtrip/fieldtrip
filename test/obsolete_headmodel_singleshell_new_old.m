function test_headmodel_singleshell_new_old

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_headmodel_singleshell ft_prepare_headmodel ft_headmodel_singleshell

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

grad.coilpos = pnt * 120;
grad.coilori = pnt;
grad.chanpos = pnt * 120;
grad.chanori = pnt;
grad.tra = eye(size(pnt,1));
for i=1:size(pnt,1)
  grad.label{i} = sprintf('%d', i);
end

arg(1).name = 'conductivity';
arg(1).value = {[], 0.10, 0.33, 1.00};

optarg = constructalloptions(arg);

for i=1:size(optarg,1)
  
  arg = optarg(i,:);
  
  % new way - low level:
  vol{1} = ft_headmodel_singleshell(geom.bnd(1),arg{:});
  
  % old way:
  tmpcfg = ft_keyval2cfg(arg{:});
  tmpcfg.headshape = geom.bnd(1);
  vol{2} = ft_prepare_singleshell(tmpcfg);
  
  % new way - high level:
  tmpcfg = ft_keyval2cfg(arg{:});
  tmpcfg.method = 'singleshell';
  vol{3} = ft_prepare_headmodel(tmpcfg,geom.bnd(1));
  
  % compute the leadfields for a comparison
  [vol{1}, grad] = ft_prepare_vol_sens(vol{1}, grad);
  [vol{2}, grad] = ft_prepare_vol_sens(vol{2}, grad);
  [vol{3}, grad] = ft_prepare_vol_sens(vol{3}, grad);
  lf{1} = ft_compute_leadfield([0 10 60], grad, vol{1});
  lf{2} = ft_compute_leadfield([0 10 60], grad, vol{2});
  lf{3} = ft_compute_leadfield([0 10 60], grad, vol{3});
  
  % compare the leadfields in all possible combinations
  comb = nchoosek(1:numel(vol),2);
  for j=1:size(comb,1)
    chk = comb(j,:);
    err = norm(lf{chk(1)} - lf{chk(2)}) / norm(lf{chk(1)});
    if err>0.001
      error('combination %d %d not successful\n',chk(1),chk(2));
    end
  end
  
end
