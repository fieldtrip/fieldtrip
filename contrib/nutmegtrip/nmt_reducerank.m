function leadgrid = nmt_reducerank(cfg,leadgrid)
% leadgrid = nmt_reducerank(cfg,leadgrid);
%
% Reduces the rank of a lead field after it has been computed. Particularly
% useful for testing reduced rank versions of BEM/FEM models.
% e.g., use cfg.reducerank = 2;


inside_idx = find(leadgrid.inside);
for ii=1:length(inside_idx)
    leadgrid.leadfield{inside_idx(ii)} = nmt_svdtrunc(leadgrid.leadfield{inside_idx(ii)},1:cfg.reducerank);
end