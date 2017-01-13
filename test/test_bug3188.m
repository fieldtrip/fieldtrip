function test_bug3188

% WALLTIME 00:30:00
% MEM 8gb

% TEST ft_connectivityanalysis

%%

global ft_default

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3188.mat'))

%%
% first issue: inside definition is still old style
cfg = [];
cfg.method = 'coh';
coh = ft_connectivityanalysis(cfg, data);
assert(all(islogical(coh.inside)));
assert(~isfield(coh, 'outside'));

% second issue: powcorr_ortho gives strangely formatted output
cfg = [];
cfg.method = 'powcorr_ortho';
pco = ft_connectivityanalysis(cfg, data);
assert(size(pco.powcorrspctrm,1)==size(pco.pos,1));


