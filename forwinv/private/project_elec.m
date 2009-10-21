function [el] = project_elec(elc, pnt, tri)

% PROJECT_ELEC projects electrodes on a triangulated surface
% and returns triangle index, la/mu parameters and distance
%
% [el] = project_elec(elc, pnt, tri)
%
% it returns a Nx4 matrix with [tri, la, mu, dist] for each electrode
%
% See also TRANSFER_ELEC 

% Copyright (C) 1999-2002, Robert Oostenveld
%
% $Log: project_elec.m,v $
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.2  2003/03/11 14:45:37  roberto
% updated help and copyrights
%

Nelc = size(elc,1);
Npnt = size(pnt,1);
Ntri = size(tri,1);
el   = zeros(Nelc, 4);

for i=1:Nelc
  smallest_dist = Inf;

  for j=1:Ntri
    [proj, dist] = ptriproj(pnt(tri(j,1),:), pnt(tri(j,2),:), pnt(tri(j,3),:), elc(i,:), 1);
    if dist<smallest_dist
      % remember the triangle index, distance and la/mu
      [la, mu] = lmoutr(pnt(tri(j,1),:), pnt(tri(j,2),:), pnt(tri(j,3),:), proj);
      smallest_dist = dist; 
      smallest_tri  = j; 
      smallest_la   = la; 
      smallest_mu   = mu; 
    end
  end

  % store the projection for this electrode
  el(i,:) = [smallest_tri smallest_la smallest_mu smallest_dist];
end

