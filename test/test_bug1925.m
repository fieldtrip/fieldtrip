function test_bug1925

% MEM 1500mb
% WALLTIME 00:10:00

% TEST surface_nesting ft_headmodel_bemcp

[ftver, ftpath] = ft_version;
cd(fullfile(ftpath, 'forward/private')); % this is where the surface_nesting function is located

[pos, tri] = icosahedron162;

bnd10.id  = 10;
bnd10.pos = pos*10;
bnd10.tri = tri;

bnd20.id  = 20;
bnd20.pos = pos*20;
bnd20.tri = tri;

bnd30.id  = 30;
bnd30.pos = pos*30;
bnd30.tri = tri;

bnd40.id  = 40;
bnd40.pos = pos*40;
bnd40.tri = tri;

bnd50.id  = 50;
bnd50.pos = pos*50;
bnd50.tri = tri;

bnd    = bnd10;
bnd(2) = bnd20;
bnd(3) = bnd30;
bnd(4) = bnd40;
bnd(5) = bnd50;
assert(equalorder(surface_nesting(bnd, 'insidefirst'), 1:5));
assert(equalorder(surface_nesting(bnd, 'outsidefirst'), fliplr(1:5)));

bnd    = bnd10;
bnd(3) = bnd20;
bnd(2) = bnd30;
bnd(5) = bnd40;
bnd(4) = bnd50;
assert(equalorder(surface_nesting(bnd, 'insidefirst'), [1 3 2 5 4]));
assert(equalorder(surface_nesting(bnd, 'outsidefirst'), fliplr([1 3 2 5 4])));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to deal with row and column comparisons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = equalorder(a, b)
bool = isequal(a(:), b(:));

