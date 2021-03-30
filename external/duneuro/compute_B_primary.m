% compute primary magnetic B-field analytically
%
% input:
% coils (Nx3 matrix)
% dipoles (Mx6 matrix)
% projections (Nx3) matrix)

function [Bp] = compute_B_primary(coils, dipoles, projections)

%check input
if size(coils,2)~=3
  error('Column size of coils must be 3.')
end

if size(dipoles,2)~=6
  error('Column size of dipoles must be 6.')
end

if size(projections,2)~=3
  error('Column size of projections must be 3.')
end

% apply formula of Sarvas

dip_pos = dipoles(:,1:3);
dip_mom = dipoles(:,4:6);
Bp = zeros(size(coils,1), size(dipoles,1));
for i = 1:size(coils,1)
  for j = 1 : size(dip_pos,1)
    R = coils(i,:);
    R_0 = dip_pos(j,:);
    A = R - R_0;
    a = norm(A);
    aa = A./(a^3);
    
    BpHelp = cross(dip_mom(j,:),aa);
    Bp(i,j) = BpHelp * projections(i, :)'; % projection of the primary B-field along the coil orientations
  end
end