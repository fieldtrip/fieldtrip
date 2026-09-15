function [Bx, By, Bz] = circularcoil(x, y, z, d, I)

% CIRCULARCOIL magnetic field of a circular current loop
%
% Use as:
%   [Bx, By, Bz] = circularcoil(x, y, z, d, I)
%
% where
%   x,y,z : observation point(s) [m]
%   d     : coil diameter [m]
%   I     : coil current [A]
% and
%   Bx,By,Bz : magnetic field components [T]
%
% The coil is centered at (0,0,0), lies in the xy-plane,
% and carries current I. Positive I gives a field in +z
% at the center.

    mu0 = 4*pi*1e-7;
    a = d/2;

    % Cylindrical radial coordinate
    rho = sqrt(x.^2 + y.^2);

    % Allocate output
    Bx = zeros(size(x));
    By = zeros(size(x));
    Bz = zeros(size(x));

    % ------------------------------------------------------------
    % Points away from the symmetry axis
    % ------------------------------------------------------------
    ind = rho > 0;

    r = rho(ind);
    zz = z(ind);

    % Elliptic integral parameter
    m = 4*a*r ./ ((a+r).^2 + zz.^2);

    % Complete elliptic integrals
    [K,E] = ellipke(m);

    Delta = sqrt((a+r).^2 + zz.^2);

    denom = (a-r).^2 + zz.^2;

    % Cylindrical components
    Brho = mu0*I*zz ./ (2*pi*r.*Delta) .* ...
        (-K + ...
        (a^2 + r.^2 + zz.^2) ./ denom .* E);

    Bzz = mu0*I ./ (2*pi*Delta) .* ...
        (K + ...
        (a^2 - r.^2 - zz.^2) ./ denom .* E);

    % Convert cylindrical -> Cartesian
    Bx(ind) = Brho .* x(ind) ./ r;
    By(ind) = Brho .* y(ind) ./ r;
    Bz(ind) = Bzz;

    % ------------------------------------------------------------
    % Points on the symmetry axis
    % ------------------------------------------------------------
    ind0 = ~ind;

    Bz(ind0) = mu0*I*a^2 ./ ...
        (2*(a^2 + z(ind0).^2).^(3/2));

end
