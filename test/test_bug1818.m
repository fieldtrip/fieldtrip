function test_bug1818

% test the reading function of meshes used for constructing SIMBIO FEM head models
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1818

% TEST test_bug1818
% TEST ft_read_headshape ft_datatype_parcellation

% vista
mesh = ft_read_headshape('/home/common/matlab/fieldtrip/data/test/bug1818/cube2mm3layervorwerk_ns_127_127_127.v');

parcellation1 = ft_datatype_parcellation(mesh);
parcellation2 = ft_datatype_parcellation(mesh,'parcellationstyle','probabilistic');

if~(ft_datatype(parcellation1,'parcellation'))
  error('the conversion to a parcellation failed');
end

if~(ft_datatype(parcellation2,'parcellation'))
  error('the conversion to a parcellation failed');
end

% tetgen-elements
mesh = ft_read_headshape('/home/common/matlab/fieldtrip/data/test/bug1818/tet_4layer_127_127_127.1.ele');

parcellation1 = ft_datatype_parcellation(mesh);
parcellation2 = ft_datatype_parcellation(mesh,'parcellationstyle','probabilistic');

if~(ft_datatype(parcellation1,'parcellation'))
  error('the conversion to a parcellation failed');
end

if~(ft_datatype(parcellation2,'parcellation'))
  error('the conversion to a parcellation failed');
end

% tetgen-node
mesh = ft_read_headshape('/home/common/matlab/fieldtrip/data/test/bug1818/tet_4layer_127_127_127.1.node');

parcellation1 = ft_datatype_parcellation(mesh);
parcellation2 = ft_datatype_parcellation(mesh,'parcellationstyle','probabilistic');

if~(ft_datatype(parcellation1,'parcellation'))
  error('vistaheadshapeisnotparcellationdatatype');
end

if~(ft_datatype(parcellation2,'parcellation'))
  error('vistaheadshapeisnotparcellationdatatype');
end


