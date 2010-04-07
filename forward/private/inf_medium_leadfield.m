function [lf] = inf_medium_leadfield(rd, pnt, cond);

% INF_MEDIUM_LEADFIELD calculate the infinite medium leadfield
% on positions pnt for dipole position R and conductivity cond
%       
% [lf] = inf_medium_leadfield(R, pnt, cond)

% Copyright (C) 1998, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

