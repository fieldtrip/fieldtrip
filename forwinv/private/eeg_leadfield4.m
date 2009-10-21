function [lf, vol] = eeg_leadfield4(R, elc, vol)

% EEG_LEADFIELD4 electric leadfield for a dipole in 4 concentric spheres
% 
% [lf] = eeg_leadfield4(R, elc, vol)
%
% with input arguments
%   R	         position of the dipole
%   elc	       position of the electrodes
% and vol being a structure with the elements
%   vol.r      radius of the 4 spheres 
%   vol.c      conductivity of the 4 spheres
%   vol.t      constant factors for series expansion (optional)
%
% See also EEG_LEADFIELD4_PREPARE for precomputing the constant factors,
% which can save time when multiple leadfield computations are done.

% Copyright (C) 2002, Robert Oostenveld
%
% this implementation is adapted from
%   Lutkenhoner, Habilschrift 1992.
% the original reference is
%  Cuffin BN, Cohen D. Comparison of the magnetoencephalogram and electroencephalogram. Electroencephalogr Clin Neurophysiol. 1979 Aug;47(2):132-46. 
%
% $Log: eeg_leadfield4.m,v $
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.8  2008/12/24 13:33:27  roboos
% changed some & and | into && and ||
%
% Revision 1.7  2006/05/01 08:13:51  roboos
% fixed literature reference
%
% Revision 1.6  2003/12/05 09:48:45  roberto
% disabled error check for dipole outside brain [for EEGLAB]
%
% Revision 1.5  2003/12/04 10:43:51  roberto
% added error when dipole is outside brain compartment
%
% Revision 1.4  2003/07/29 16:04:55  roberto
% fixed bug in determining whether electrodes are lying on sphere surface
%
% Revision 1.3  2003/07/29 15:52:44  roberto
% fixed a bug in the implementation of eeg_leadfield4, caused by mixing the constants of
% lutkenhoner and cuffin
% furthermore multiple cosmetic changes and default projection of electrodes to sphere
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

% sort the spheres from the smallest to the largest
[vol.r, indx] = sort(vol.r);
[vol.c]       = vol.c(indx);

% use more convenient names for the radii and conductivity
r1 = vol.r(1); c1 = vol.c(1);
r2 = vol.r(2); c2 = vol.c(2);
r3 = vol.r(3); c3 = vol.c(3);
r4 = vol.r(4); c4 = vol.c(4);

% check whether the electrode ly on the sphere, allowing 0.5% tolerance
dist = sqrt(sum(elc.^2,2));
if any(abs(dist-r4)>r4*0.005)
  warning('electrodes do not ly on sphere surface -> using projection')
end
elc = r4 * elc ./ [dist dist dist];

% check whether the dipole is inside the brain [disabled for EEGLAB]
% if sqrt(sum(R.^2))>=r1
%   error('dipole is outside the brain compartment');
% end

% rotate everything so that the dipole is along the pos. z-axis
% only if the dipole is not in the origin or along the positive z-axis
if R(1)~=0 | R(2)~=0
  % compute the rotation matrix
  % the inverse rotation matrix is the transposed of this one
  val1 = norm(R);
  val2 = norm(R(1:2));
  rot(1,1) = R(1) * R(3) / (val1 * val2); 
  rot(1,2) = R(2) * R(3) / (val1 * val2);
  rot(1,3) = -1.0 * val2 / val1;
  rot(2,1) = -1.0 * R(2) / val2;
  rot(2,2) =        R(1) / val2;
  rot(2,3) =                  0; 
  rot(3,:) = R ./ val1;
  % rotate the electrodes
  elc = elc*rot';
elseif R(1)==0 && R(2)==0 && R(3)<0
  % dipole on negative z-axis, rotation is very simple: around x-axis
  elc(2,:) = -elc(2,:);
  elc(3,:) = -elc(3,:);
else
  % dipole is on positive z-axis, nothing has to be done
end

% compute the constant factors for the sphere configuration if needed
if ~isfield(vol, 't')
  vol.t = eeg_leadfield4_prepare(vol);
end

Nchans = size(elc,1);
lf     = zeros(Nchans,3);
Nmax   = length(vol.t);
n      = 1:Nmax;
f      = norm(R)/r4;		% following cuffin1979
% c      = r2/r4;		% following cuffin1979
% d      = r3/r4;		% following cuffin1979

% this code is to cross-validate the lutkenhoner and cuffin implementations
% [lut_t, cuf_t] = eeg_leadfield4_prepare(vol);
% lut_c = (2*n+1).^4.*f.^(n-1) ./ (lut_t.*4*pi*c4*r4^2);
% cuf_c = (2*n+1).^4.*f.^(n-1) .*(c*d).^(2.*n+1) ./ (cuf_t.*4*pi*c4*r4^2);

% given a fixed volume conductor, these only need to be computed once for all electrodes
const = (2*n+1).^4.*f.^(n-1) ./ (vol.t.*4*pi*c4*r4^2);

for i=1:Nchans
  % convert the position of the electrodes to spherical coordinates
  [phi, el] = cart2sph(elc(i,1), elc(i,2), elc(i,3));

  % change from colatitude to latitude and compute the cosine 
  cos_theta = cos(pi/2-el);

  % the series summation starts at zero
  s_x = 0;
  s_z = 0;

  for n=1:Nmax
    P0  = plgndr(n,0,cos_theta);		% zero'th order Legendre
    P1  = plgndr(n,1,cos_theta);		% first order Legendre
    s_x = s_x + const(n)*P1/n;			% s_y is identical
    s_z = s_z + const(n)*P0;
  end

  lf(i,1) = -cos(phi) * s_x;
  lf(i,2) = -sin(phi) * s_x;			% s_y is identical to s_x
  lf(i,3) = 1         * s_z;
end

% apply the inverse rotation to the leadfield matrix
if R(1)~=0 || R(2)~=0
  lf = lf*rot;
elseif R(1)==0 && R(2)==0 && R(3)<0
  lf(2,:) = -lf(2,:);
  lf(3,:) = -lf(3,:);
end

