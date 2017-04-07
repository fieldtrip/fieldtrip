function b=om_spher_mag_iso(x,q,p)
% Calculates the magnetic field (a vector) at point x on the surface of a layered
% conductive sphere induced by a current dipole of strength q at point
% p. The norm of x must be equal to the last rk, at least.
%
% All vectors (x,q,p) should be in 3D.
%
% Based on
% MOSHER & Al. : EEG and MEG: Forward solutions for inverse methods
% Ilmoniemi & Al
% Sarvas
% 
% Author: E. Olivi 2008/11/18

% Copyright (C) 2010-2017, OpenMEEG developers

mu0 = 4*pi*1e-7;
%mu0 = 1;

[n,m] = size(x);
if m~=3,
  error('x should be a set of row 3D vectors.');
end

b=zeros(n,3);
for ipos=1:n;
    xi=x(ipos,:);
    xx=norm(xi);
    d=xi-p;
    dd=norm(d);

    F=dd*(xx*dd+xx^2-dot(xi,p));

    gradF=(dd^2/xx+dot(d,xi)/dd+2*dd+2*xx)*xi-(dd+2*xx+dot(d,xi)/dd)*p;

    b(ipos,:)=mu0/(4*pi*F^2)*(F*cross(q,p)-(dot(cross(q,p),xi)*gradF));
end;
