function [lf] = inf_medium_leadfield(rd, pnt, cond);

% INF_MEDIUM_LEADFIELD calculate the infinite medium leadfield
% on positions pnt for dipole position R and conductivity cond
%       
% [lf] = inf_medium_leadfield(R, pnt, cond)

% Copyright (C) 1998, Robert Oostenveld
%
% $Log: inf_medium_leadfield.m,v $
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.5  2005/02/23 14:31:19  roboos
% changed the detection of Ndipoles, added reshaping of Nx3 input
%
% Revision 1.4  2003/08/04 09:12:32  roberto
% added check for dipole on BEM model boundary, gives warning message
%
% Revision 1.3  2003/06/03 08:29:26  roberto
% fixed error affecting multiple dipole computations
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

siz = size(rd);
if any(siz==1)
  % positions are specified as a single vector
  Ndipoles = prod(siz)/3;
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Ndipoles = siz(1);
  rd = rd';
  rd = rd(:);
else
  error('incorrect specification of dipole locations');
end

Npnt     = size(pnt,1);
lf       = zeros(Npnt,3*Ndipoles);
s1       = size(rd);

if s1(1)>s1(2)
  % make sure that the dipole position is a row-vector
  rd = rd';
end

for i=1:Ndipoles
  r = pnt - ones(Npnt,1) * rd((1:3) + 3*(i-1));
  R = (4*pi*cond) * (sum(r' .^2 ) .^ 1.5)';
  if any(R)==0
    warning('dipole lies on boundary of volume model');
  end
  lf(:,(1:3) + 3*(i-1)) = r ./ [R R R];
end

