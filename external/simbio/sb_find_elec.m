function diri = sb_find_elec(vol,elc)
diri = zeros(size(elc,1),1);
for i=1:size(elc,1)
    [dist, diri(i)] = min(sum(bsxfun(@minus,vol.wf.nd,elc(i,:)).^2,2));
end
end