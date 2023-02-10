function test_bug1992

% DEPENDENCY ft_filetype ft_read_headshape ft_write_headshape loadjson savejson loadbj savebj jsonopt

filename = [tempname '.jmsh'];

[pos, tri] = mesh_sphere(162);

% combine them in a structure
mesh1.pos = pos;
mesh1.tri = int32(tri);
mesh1.info = struct('JMeshVersion', '0.5', 'Dimension', 3, ...
                    'AnnotationFormat', 'https://github.com/NeuroJSON/jmesh/blob/master/JMesh_specification.md', ...
                    'SerialFormat', 'http://json.org');
mesh1.poslabel = ones(1, size(pos,1));
mesh1.trilabel = ones(1, size(tri,1));

ft_write_headshape(filename, mesh1, 'format', 'neurojson_jmesh');
mesh2 = ft_read_headshape(filename);
assert(isequal(mesh1.pos, mesh2.pos));
assert(isequal(mesh1.poslabel, mesh2.poslabel));
assert(isequal(mesh1.tri, mesh2.tri));
assert(isequal(mesh1.trilabel, mesh2.trilabel));
delete(filename);

ft_write_headshape(filename, mesh1, 'format', 'neurojson_jmesh', 'jsonopt', {'compact',1, 'compression', ''});
mesh2 = ft_read_headshape(filename);
assert(isequal(mesh1.tri, mesh2.tri));
ratio = mesh1.pos(:)\mesh2.pos(:);
assert(ratio>0.99 && ratio<1.01);
assert(isequal(mesh1.poslabel, mesh2.poslabel));
assert(isequal(mesh1.trilabel, mesh2.trilabel));
delete(filename);

filename = [tempname '.bmsh'];
mesh1.info.SerialFormat = 'https://neurojson.org/bjdata/draft2';

ft_write_headshape(filename, mesh1, 'format', 'neurojson_bmesh');
mesh2 = ft_read_headshape(filename);
assert(isequal(mesh1.pos, mesh2.pos));
assert(isequal(mesh1.poslabel, mesh2.poslabel));
assert(isequal(mesh1.tri, mesh2.tri));
assert(isequal(mesh1.trilabel, mesh2.trilabel));
delete(filename);





