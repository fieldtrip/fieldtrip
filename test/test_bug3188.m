function test_bug3188

% WALLTIME 00:30:00
% MEM 6gb
% DEPENDENCY ft_connectivityanalysis
% DATA private

%%

load(dccnpath('/project/3031000.02/test/bug3188.mat'))

%%
% first issue: inside definition is still old style

cfg = [];
cfg.method = 'coh';
coh = ft_connectivityanalysis(cfg, data);
assert(all(islogical(coh.inside)));
assert(~isfield(coh, 'outside'));

%%
% second issue: powcorr_ortho gives strangely formatted output

cfg = [];
cfg.method = 'powcorr_ortho';
pco = ft_connectivityanalysis(cfg, data);
assert(size(pco.powcorrspctrm,1)==size(pco.pos,1));


