function source = nmt_polaritytweak(cfg,source)

inside_idx = find(source.inside);

s = cell2mat(source.avg.mom(inside_idx));
ori = cell2mat(source.avg.ori(inside_idx));

if(isfield(cfg,'toilim'))
    tsel = dsearchn(source.time',cfg.toilim');
    tsel = tsel(1):tsel(2);
else
    [dum,tsel]=max(max(abs(s)));
end

% flip based on polarity of voxel with maximum power in desired time window
[dum,b]=max(nmt_rownorm(s(:,tsel)));
flipper=sign(s(b,tsel)*s(:,tsel)').*speye(size(s,1)); 
s = flipper*s;
ori = ori*flipper; 

for ii=1:size(s,1)
    source.avg.mom{inside_idx(ii)} = s(ii,:);
    source.avg.ori{inside_idx(ii)} = ori(:,ii);
end
