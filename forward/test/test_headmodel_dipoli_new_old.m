function test_headmodel_dipoli_new_old

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

arg(1).name = 'conductivity';
arg(2).name = 'isolatedsource';

arg(1).value = {[], [1 1/20 1], [0.33 0.125 0.33], [1 1 1], [0.1 0.1 0.1]};% [], 
arg(2).value = {'yes' , 'no'};

optarg = constructalloptions(arg);

for i=1:size(optarg,1)
  
  arg = optarg(i,:);
  
  % new way - low level:
  vol1 = ft_headmodel_bem_dipoli(geom,arg{:});

  % old way:
  tmpcfg = keyval2cfg(arg{:});
  tmpcfg.method = 'dipoli';
  vol2 = ft_prepare_bemmodel(tmpcfg,geom);

  % new way - high level:
  tmpcfg = keyval2cfg(arg{:});
  tmpcfg.method = 'bem_dipoli';
  vol3 = ft_prepare_headmodel(tmpcfg,geom.bnd); % FIXME: this argument has to be the same

  % compute an example leadfield
  if ~isequal(vol1,vol2) || ~isequal(vol2,vol3) || ~isequal(vol1,vol3)
    error('not successful')
  end
  
end

return

% create a set of electrodes
sel = find(pnt(:,3)>0);
sens.pnt = pnt(sel,:) * 100;
for i=1:length(sel)
  sens.label{i} = sprintf('chan%03d', i);
end

% project the electrodes on the volume conduction model
[vol1, sens] = ft_prepare_vol_sens(vol1, sens);
[vol2, sens] = ft_prepare_vol_sens(vol2, sens);

dip = {
  [ 0 0 50]
  };


for i=1:numel(dip)
  % conpute an example leadfield
  lf1 = ft_compute_leadfield(dip{i}, sens, vol1);
  lf2 = ft_compute_leadfield(dip{i}, sens, vol2);
  
  if norm(lf1-lf2)/(norm(lf1+lf2)/2)>2*eps
    error('not successful')
  end
end




