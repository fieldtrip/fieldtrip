%  HBF_CHECKMESH performs some tests on a boundary mesh
% 
%  status=HBF_CHECKMESH(mesh,lasttesttodo)
%  status=HBF_CHECKMESH(points,elements,lasttesttodo)
%    mesh:   hbf struct for triangle mesh
%    points: mesh vertices, [N x 3]
%    elements:   mesh triangle description, [M x 3]
%    lasttesttodo: last test to do --- between 1 and 8.
% 
%  Tests 1--4 are basic tests that should be carried out always.
%  Tests 5 -- 7 test for a closed mesh. If some of tests 2--4 fail, 5--7
%  will probably not work correctly.
% 
%  test 1: does each vertex belong to some triangle?
%  test 2: does each triangle consist of three different vertices?
%  test 3: does each triplet of vertices form one and only one triangle?
%  test 4: does each vertex have a unique position?
%  test 5: does each vertex belong to at least 3 triangles?
%  test 6: is the solid angle spanned by the mesh at a point outside the mesh 0?
%  test 7: does each triangle-side belong to exactly two triangles?
% 
%  This tool only spots basic errors; it does not attempt to correct them.
% 
%  v180615 (c) Matti Stenroos
%