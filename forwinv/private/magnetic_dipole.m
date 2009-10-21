function [lf] = magnetic_dipole(R, pos, ori)

% MAGNETIC_DIPOLE leadfield for a magnetic dipole in an infinite medium
%
% [lf] = magnetic_dipole(R, pos, ori)
%
% with input arguments
%   R           position dipole
%   pos         position magnetometers
%   ori         orientation magnetometers

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: magnetic_dipole.m,v $
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.1  2003/03/12 09:22:18  roberto
% new implementation, optimized for speed
% based on http://scienceworld.wolfram.com/physics/MagneticDipole.html
%

u0 = 1e-7;
nchan = size(pos,1);

% ensure that the dipole position is a row vector
R = reshape(R, [1 3]);

% shift the magnetometers so that the dipole is in the origin
pos = pos - repmat(R, [nchan 1]);

% change the variable names for convenience
% R position of magnetometer, relative to dipole
% r distance of magnetometer, relative to dipole
R = pos;
r = sqrt(sum(R.^2,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the original code, which follows the physical formulation closely
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnetic field on all magnetometer positions for an x-oriented dipole
% Bmx = u0/(4*pi) * (3 * repmat(R(:,1), [1 3]) .* R - repmat([1 0 0], [nchan 1]) .* repmat(r.^2, [1 3])) ./ repmat(r.^5, [1 3]);
% magnetic field on all magnetometer positions for an y-oriented dipole
% Bmy = u0/(4*pi) * (3 * repmat(R(:,2), [1 3]) .* R - repmat([0 1 0], [nchan 1]) .* repmat(r.^2, [1 3])) ./ repmat(r.^5, [1 3]);
% magnetic field on all magnetometer positions for an z-oriented dipole
% Bmz = u0/(4*pi) * (3 * repmat(R(:,3), [1 3]) .* R - repmat([0 0 1], [nchan 1]) .* repmat(r.^2, [1 3])) ./ repmat(r.^5, [1 3]);
% compute the field along the orientation of each magnetometer for an x/y/z oriented dipole
% lf(:,1) = dot(Bmx, ori, 2);
% lf(:,2) = dot(Bmy, ori, 2);
% lf(:,3) = dot(Bmz, ori, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the optimized code, which make the previous computation approximately
% three times more efficent by introducing temporary variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r2 = repmat(r.^2, [1 3]);
r5 = repmat(r.^5, [1 3]);
x = R(:,1); x = [x x x];
y = R(:,2); y = [y y y];
z = R(:,3); z = [z z z];
mx = zeros(nchan,3); mx(:,1) = 1;
my = zeros(nchan,3); my(:,2) = 1;
mz = zeros(nchan,3); mz(:,3) = 1;
Tx = (3 * x .* R - mx .* r2);
Ty = (3 * y .* R - my .* r2);
Tz = (3 * z .* R - mz .* r2);
lf(:,1) = dot(Tx, ori, 2);
lf(:,2) = dot(Ty, ori, 2);
lf(:,3) = dot(Tz, ori, 2);
lf = u0/(4*pi) * lf ./ r5;

