function failed_bug2790

% MEM 6gb
% WALLTIME 00:10:00

% TEST test_bug2790
% TEST ft_selectdata
% TEST ft_connectivityanalysis

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2790.mat');
load(filename);

% something seems to be strange with the data, don't worry about that for
% now
in  = find(data.inside);
sel = ~cellfun('isempty',data.avg.mom(in));
in  = in(sel);
data.inside = false(size(data.pos,1),1);
data.inside(in) = true;
if isfield(data, 'outside'), data = rmfield(data, 'outside'); end

rpt = cfg.trials;


cfg.trials = 'all';
coh_all = ft_connectivityanalysis(cfg,data);

mom    = cat(1,data.avg.mom{in});
refmom = mom(in==cfg.refindx,:);
coh    = abs(imag(mom*refmom'))./sqrt(sum(abs(mom).^2,2).*sum(abs(refmom).^2,2));

assert(norm(coh_all.cohspctrm(in)-coh)<1e-15);

cfg.trials = rpt;
coh_rpt = ft_connectivityanalysis(cfg,data);
