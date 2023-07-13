function test_bug1992

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_filetype ft_read_headshape ft_write_headshape loadjson savejson loadbj savebj jsonopt
% DATA private

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/mesh/jmesh'));

%%
% see https://github.com/NeuroJSON/JMeshSamples

filename = {
  'surface/happy_budda.jmsh'
  % 'surface/happy_budda_lzma.jmsh'
  'surface/happy_budda_zlib.jmsh'
  'surface/porous_tri.jmsh'
  'surface/sidecut_fiber_tri.jmsh'
  'surface/skull_tri_multipart_by_name_zlib.jmsh'
  'surface/stanford_bunny.jmsh'
  'test/cube_quad.jmsh'
  'test/cube_tri.jmsh'
  'test/cube_tri_annotated_array.jmsh'
  'test/cube_tri_zlib.jmsh'
  'test/cyl_plc.jmsh'
  'test/isosphere_tet.jmsh'
  'test/isosphere_tri.jmsh'
  'test/mobius_quad.jmsh'
  'test/mobius_tri.jmsh'
  'test/sidecut_fiber_plc.jmsh'
  'test/sphere_quad.jmsh'
  'test/sphere_tri.jmsh'
  % 'test/twocube_csg_union.jmsh'
  'test/twocube_plc.jmsh'
  'tetmesh/dumbbell.jmsh'
  'tetmesh/sidecut_fiber_tet_flex.jmsh'
  'tetmesh/skull_tet.jmsh'
  'tetmesh/sphbox_tet_flex.jmsh'
  'tetmesh/stanford_bunny_tet.jmsh'
  };


outputfile = tempname;

for i=1:numel(filename)
  disp('----------------')
  disp(filename{i})
  mesh = ft_read_headshape(filename{i});

  figure
  ft_plot_mesh(mesh)
  camlight

  ft_write_headshape([outputfile '.jmsh'], mesh);
  delete([outputfile '.jmsh']);

  ft_write_headshape([outputfile '.bmsh'], mesh);
  delete([outputfile '.msh']);
end


%%

[pos, tri] = mesh_sphere(162);

% combine them in a structure
mesh1 = [];
mesh1.info = struct('JMeshVersion', '0.5', 'Dimension', 3, ...
  'AnnotationFormat', 'https://github.com/NeuroJSON/jmesh/blob/master/JMesh_specification.md', ...
  'SerialFormat', 'http://json.org');

mesh1.pos  = pos;
mesh1.tri  = tri;
mesh1.line = tri(1:2:end,[1 2]); % keep one edge of each second triangle

% probabilistic parcellations of the vertices, triangles and lines
mesh1.Brodmann_Area_1 = rand(size(mesh1.pos,1), 1);
mesh1.Brodmann_Area_2 = rand(size(mesh1.tri,1), 1);
mesh1.Brodmann_Area_3 = rand(size(mesh1.line,1), 1);

%%

figure; ft_plot_mesh(keepfields(mesh1, {'pos'}));

%%

figure; ft_plot_mesh(keepfields(mesh1, {'pos', 'tri'}), 'vertexcolor', mesh1.Brodmann_Area_1);
figure; ft_plot_mesh(keepfields(mesh1, {'pos', 'tri'}), 'facecolor', mesh1.Brodmann_Area_2);

%%

figure; ft_plot_mesh(keepfields(mesh1, {'pos', 'line'}));

%%

filename = [tempname '.jmsh'];

ft_write_headshape(filename, mesh1, 'format', 'neurojson_jmesh');
mesh2 = ft_read_headshape(filename);
assert(isequal(mesh1.pos, mesh2.pos));
assert(isequal(mesh1.tri, mesh2.tri));
assert(isequal(mesh1.line, mesh2.line));
delete(filename);

%%

ft_write_headshape(filename, mesh1, 'format', 'neurojson_jmesh', 'jmeshopt', {'compact',1, 'compression', ''});
mesh2 = ft_read_headshape(filename, 'jmeshopt', {'formatversion', 2});
assert(isequal(mesh1.tri, mesh2.tri));
ratio = mesh1.pos(:)\mesh2.pos(:);
assert(ratio>0.99 && ratio<1.01);
assert(isequal(mesh1.line, mesh2.line));
delete(filename);

%%

filename = [tempname '.bmsh'];
mesh1.info.SerialFormat = 'https://neurojson.org/bjdata/draft2';

ft_write_headshape(filename, mesh1, 'format', 'neurojson_bmesh');
mesh2 = ft_read_headshape(filename);
assert(isequal(mesh1.pos, mesh2.pos));
%assert(isequal(mesh1.poslabel, mesh2.poslabel));
assert(isequal(mesh1.tri, mesh2.tri));
%assert(isequal(mesh1.trilabel, mesh2.trilabel));
assert(isequal(mesh1.line, mesh2.line));
%assert(isequal(mesh1.poly, mesh2.poly));
delete(filename);
