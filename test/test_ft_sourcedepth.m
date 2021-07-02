function test_ft_sourcedepth

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourcedepth

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

% A negative depth indicates that the source is inside the source
% compartment, positive indicates outside.

assert(isalmostequal(ft_sourcedepth([  0   0  40], headmodel), -120));
assert(isalmostequal(ft_sourcedepth([  0   0  50], headmodel), -110));
assert(isalmostequal(ft_sourcedepth([  0   0  60], headmodel), -100));
assert(isalmostequal(ft_sourcedepth([  0   0 160], headmodel),    0));
assert(isalmostequal(ft_sourcedepth([  0   0 200], headmodel),   40));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testMesh(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pos, tri] = mesh_sphere(2560);

headmodel = [];
headmodel.bnd.pos = pos*120;
headmodel.bnd.pos(:,3) = headmodel.bnd.pos(:,3) + 40;
headmodel.bnd.tri = tri;
headmodel.unit = 'mm';

% the volume conduction model is not a perfect sphere, hence some tolerance is needed
assert(isalmostequal(ft_sourcedepth([  0   0  40], headmodel), -120, 'abstol', 1));
assert(isalmostequal(ft_sourcedepth([  0   0  50], headmodel), -110, 'abstol', 1));
assert(isalmostequal(ft_sourcedepth([  0   0  60], headmodel), -100, 'abstol', 1));
assert(isalmostequal(ft_sourcedepth([  0   0 160], headmodel),    0, 'abstol', 1));
assert(isalmostequal(ft_sourcedepth([  0   0 200], headmodel),   40, 'abstol', 1));
