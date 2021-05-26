function test_issue1711

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_source2full ft_source2sparse

%%

% R = rotate([10, 20, 30])
% T = translate([10, 20, 30])

% the rotation should not be so large to tilt the axes more than 45 degrees
% since the reconstructed positions will be aligned along x, y, z
R = [
   0.8138   -0.4698    0.3420         0
   0.5438    0.8232   -0.1632         0
  -0.2049    0.3188    0.9254         0
   0         0         0              1.0000
  ];


T = [
  1     0     0    10
  0     1     0    20
  0     0     1    30
  0     0     0     1
  ];

dim = [5, 5, 5];
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
transform = T*R;
pos = ft_warp_apply(transform, [X(:) Y(:) Z(:)]);

inside = true(dim);

% flag the corner points as outside
inside(1,1,1) = false;
inside(dim(1),1,1) = false;
inside(1,dim(2),1) = false;
inside(dim(1),dim(2),1) = false;
inside(1,1,dim(3)) = false;
inside(dim(1),1,dim(3)) = false;
inside(1,dim(2),dim(3)) = false;
inside(dim(1),dim(2),dim(3)) = false;

%%
% this is a volumetric description, which means it must be simple and restricted to 3D arrays

source1 = [];
source1.dim = dim;
source1.transform = transform;
source1.avg.pow = randn(dim);

% this adds the inside field
source1 = ft_checkdata(source1, 'insidestyle', 'logical');
source1.inside(~inside) = false;

source1s = ft_source2sparse(source1);
source1f = ft_source2full(source1s);

difference = round(source1f.pos - pos, 3, 'significant');
assert(max(difference(:))<1e-6);

assert(all(source1s.inside));
assert(isequal(source1f.inside, inside(:)));

assert(~isfield(source1s, 'outside'));
assert(~isfield(source1f, 'outside'));

%%
% this is a source description, which means it can be more complex

source2 = [];
source2.pos = pos;
source2.avg.pow = randn(prod(dim),1);
source2.cohspctrm = randn(prod(dim));

source2.time = 1:20;

source2.mom = repmat({randn(3,20)}, prod(dim), prod(dim), 1);
source2.mom(~inside, ~inside) = {[]};

% this adds the inside field, assuming that all sources are inside
source2 = ft_checkdata(source2, 'insidestyle', 'logical');
source2.inside(~inside) = false;

source2s = ft_source2sparse(source2);
source2f = ft_source2full(source2s);

difference = round(source2f.pos - pos, 3, 'significant');
assert(max(difference(:))<1e-6);

assert(all(source1s.inside));
assert(isequal(source1f.inside, inside(:)));

assert(~isfield(source1s, 'outside'));
assert(~isfield(source1f, 'outside'));

