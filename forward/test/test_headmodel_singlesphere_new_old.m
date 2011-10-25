function test_headmodel_concentricspheres_new_old

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
arg(1).value = {[], [1 1/20 1], [0.33 0.125 0.33], [1 1 1], [0.1 0.1 0.1]};

optarg = constructalloptions(arg);
% random shuffle the configurations
optarg = optarg(randperm(size(optarg,1)), :);


for i=1:size(optarg,1)
  
  arg = optarg(i,:);
  
  % new way - low level: singlesphere
  vol{1} = ft_headmodel_singlesphere(geom.bnd(1),arg{:});
  vol{1} = rmfield(vol{1},'unit');
  vol{1} = rmfield(vol{1},'type');

  % old way - low level: concentricspheres
  tmpcfg = keyval2cfg(arg{:});
  tmpcfg.headshape = geom.bnd(1);
  vol{2} = ft_prepare_concentricspheres(tmpcfg);
  vol{2} = rmfield(vol{2},'unit');
  vol{2} = rmfield(vol{2},'type');
  
  % new way - high level: singlesphere
  tmpcfg = keyval2cfg(arg{:});
  tmpcfg.method = 'singlesphere';
  vol{3} = ft_prepare_headmodel(tmpcfg,geom.bnd(1));
  vol{3} = rmfield(vol{3},'unit');
  vol{3} = rmfield(vol{3},'type');
  
  % new way - high level: concentricspheres:
  tmpcfg = keyval2cfg(arg{:});
  tmpcfg.method = 'concentricspheres';
  vol{4} = ft_prepare_headmodel(tmpcfg,geom.bnd(1));
  vol{4} = rmfield(vol{4},'type');  
  
  % compare the volume conductor structures
  comb = nchoosek(1:numel(vol),2);
  
  for i=1:size(comb,1)
    chk = comb(i,:);
    try
      if ~isequal(vol{chk(1)},vol{chk(2)})
        str = sprintf('combination %d %d not successful\n',chk(1),chk(2));
        error(str)
      end
    catch me
      fprintf(me.message)
    end
  end

end
