function diri = sb_find_elec(vol,sens)

% SB_FIND_ELEC
%
% $Id$

diri = zeros(size(sens.elecpos,1),1);
dist = zeros(size(sens.elecpos,1),1);
for i=1:size(sens.elecpos,1)
    [dist(i), diri(i)] = min(sum(bsxfun(@minus,vol.pos,sens.elecpos(i,:)).^2,2));
end
if ~all(dist < 1e-8)
    error('Electrode positions are not located on mesh nodes! This should not happen, please contact support.');
end
end
