function test_ft_appendsource

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_appendsource

gridsize = 1320;
nsamples = 20;
source = [];
source.pos = randn(gridsize,3);
source.time = 1:nsamples;
source.pow = randn(gridsize,nsamples);
source.inside = 1:gridsize/2;

% test append on rpt (only case that works yet)
source2 = source;
source2.pow = randn(gridsize,nsamples);
cfg = [];
cfg.parameter = 'pow';
cfg.appenddim = 'rpt';
sourceout = ft_appendsource(cfg, source, source2);

%return; %following cases do not work yet (cf. issue #1833)
%%
% test case of different positions, same time frames
source2 = source;
source2.pow = randn(gridsize,nsamples);
source2.pos = randn(gridsize,3);
cfg = [];
cfg.parameter = 'pow';
sourceout = ft_appendsource(cfg, source, source2);

% test case of different time frames, same positions
source2 = source;
source2.pow = randn(gridsize,nsamples);
source2.time = nsamples+1:nsamples*2;
cfg = [];
cfg.parameter = 'pow';
sourceout = ft_appendsource(cfg, source, source2);

% this should reorder the inputs
sourceout2 = ft_appendsource(cfg, source2, source);
assert(isequal(sourceout.time, sourceout2.time));

% this should throw an error
source2.time = source.time + 2;
try
  sourceout    = ft_appendsource(cfg, source, source2);
  ok = true;
catch
  ok = false;
end
assert(~ok); clear ok

% this should also throw an error
source2.time = 2.*source.time -nsamples - 0.5;
try
  sourceout    = ft_appendsource(cfg, source, source2);
  ok = true;
catch
  ok = false;
end
assert(~ok); clear ok

source = [];
[ftver, ftdir] = ft_version;
pwdir = pwd;
cd(fullfile(ftdir, 'private'));
[source.pos, source.tri] = mesh_sphere(200);
cd(pwdir);

source.pow = randn(size(source.pos,1),10);
source.freq = 1:10;

source2 = source;
source2.pos(:,3) = source.pos(:,3) + 2;

% concatenate across pos with one source containing a mesh
cfg = [];
cfg.parameter = 'pow';
sourceout = ft_appendsource(cfg, rmfield(source, 'tri'), source2);
assert(min(sourceout.tri(:))==size(source.pos,1)+1);

% concatenate across pos with 2 meshes
sourceout = ft_appendsource(cfg, source, source2);
assert(max(sourceout.tri(:))==size(sourceout.pos,1));


