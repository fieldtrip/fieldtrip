function test_ft_inside_headmodel

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_inside_headmodel

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  func = localfunctions;
  for i=1:numel(func)
    fprintf('evaluating %s\n', func2str(func{i}));
    feval(func{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testSphere(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

headmodel = [];
headmodel.r = 120;
headmodel.o = [0 0 40];
headmodel.unit = 'mm';

assert( ft_inside_headmodel([  0   0  40], headmodel));
assert( ft_inside_headmodel([119   0  40], headmodel));
assert(~ft_inside_headmodel([121   0  40], headmodel)); % outside
assert( ft_inside_headmodel([  0   0 -79], headmodel));
assert(~ft_inside_headmodel([  0   0 -81], headmodel)); % outside

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testMesh(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pos, tri] = mesh_sphere(162);

headmodel = [];
headmodel.bnd.pos = pos*120;
headmodel.bnd.pos(:,3) = headmodel.bnd.pos(:,3) + 40;
headmodel.bnd.tri = tri;
headmodel.unit = 'mm';

assert( ft_inside_headmodel([  0   0  40], headmodel));
assert( ft_inside_headmodel([119   0  40], headmodel));
assert(~ft_inside_headmodel([121   0  40], headmodel)); % outside
assert( ft_inside_headmodel([  0   0 -79], headmodel));
assert(~ft_inside_headmodel([  0   0 -81], headmodel)); % outside
