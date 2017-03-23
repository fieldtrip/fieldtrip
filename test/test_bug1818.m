function test_bug1818

% MEM 2gb
% WALLTIME 00:10:00

% test the reading function of meshes used for constructing SIMBIO FEM head models
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1818

% TEST ft_read_headshape ft_datatype_parcellation

filename = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug1818/cube2mm3layervorwerk_ns_127_127_127.v') % vista
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug1818/tet_4layer_127_127_127.1.ele')          % tetgen
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug1818/tet_4layer_127_127_127.1.node')         % tetgen
};

for i=1:length(filename)

  disp(filename{i});
  mesh = ft_read_headshape(filename{i});

  % a mesh has a pos, and possibly a tri/tet/hex field
  assert(isfield(mesh, 'pos'), 'the mesh has no pos');

  parcellation1 = ft_datatype_parcellation(mesh);
  parcellation2 = ft_datatype_parcellation(mesh,'parcellationstyle','probabilistic');

  if ~(ft_datatype(parcellation1,'parcellation'))
    error('the conversion to a parcellation failed');
  end

  if ~(ft_datatype(parcellation2,'parcellation'))
    error('the conversion to a parcellation failed');
  end

end % for each file


