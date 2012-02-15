function test_headmodel_bemcp_new_old

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

arg(1).value = {[], [1 1/20 1], [0.33 0.125 0.33], [1 1 1], [0.1 0.1 0.1]};
arg(2).value = {'yes' , 'no'};

optarg = constructalloptions(arg);
% random shuffle the configurations
optarg = optarg(randperm(size(optarg,1)), :);

for i=1:size(optarg,1)
  
  arg = optarg(i,:);
  
  % new way - low level:
  vol{1} = ft_headmodel_bemcp(geom.bnd,arg{:});

  % old way:
  tmpcfg = keyval2cfg(arg{:});
  tmpcfg.method = 'bemcp';
  vol{2} = ft_prepare_bemmodel(tmpcfg,geom);
  vol{2} = rmfield(vol{2},'unit');
  
  % new way - high level:
  tmpcfg = keyval2cfg(arg{:});
  tmpcfg.method = 'bem_cp';
  vol{3} = ft_prepare_headmodel(tmpcfg,geom.bnd);
  vol{3} = rmfield(vol{3},'unit');

  % compare the volume conductor structures
  comb = nchoosek(1:numel(vol),2);
  
  for j=1:size(comb,1)
    chk = comb(j,:);
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
