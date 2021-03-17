function test_issue1711

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_source2full ft_source2sparse

dim = [10, 11, 12];
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos = [X(:) Y(:) Z(:)];

% about 30% of the grid points will be inside
inside = rand(dim)>0.3;
outside = ~inside;

%%

source1 = [];
source1.dim = dim;
source1.transform = eye(4);
source1.pow = randn(dim);

% this adds the inside field
source1 = ft_checkdata(source1)
source1.inside(outside) = false;

source1s = ft_source2sparse(source1);
source1f = ft_source2full(source1s);

%%

source2 = [];
source2.pos = pos;
source2.dim = dim; % since it can be reshaped in a regular 3D grid
source2.pow = randn(dim);

% this adds the inside field
source2 = ft_checkdata(source2)
source2.inside(outside) = false;

source2s = ft_source2sparse(source2);
source2f = ft_source2full(source2s);

