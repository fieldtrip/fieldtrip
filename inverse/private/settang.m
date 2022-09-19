function [MDip] = settang(theta, tanu, tanv)

% set the dipole cartesian direction, given:
% 1) the instantenious decomposition vectors tanu and tanv
% 2) the instanteneous dipole orientation theta

u=cos(theta);
v=sin(theta);
for ix=1:3
  MDip(ix)=u*tanu(ix)+v*tanv(ix);
end
