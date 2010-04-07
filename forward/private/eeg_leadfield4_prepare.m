function [lut_t, cuf_t] = eeg_leadfield4_prepare(vol, Nmax);

% EEG_LEADFIELD4_PREPARE computes constant factors for series expansion
% for the 4 concentric sphere electric leadfield computation 
%
% use this function prior to repeated calls of eeg_leadfield4 according to
%   vol.t = eeg_leadfield4_prepare(vol, N);
% where
%   vol.r      radius of the 4 spheres 
%   vol.c      conductivity of the 4 spheres
% and N is the number of terms for the series (default 60). The constant
% factors t then do not have to be computed each time in eeg_leadfield4.
%
% See also EEG_LEADFIELD4

% Copyright (C) 2002, Robert Oostenveld
%
% this implementation is adapted from
%   Lutkenhoner, Habilschrift 1992.
% which again is taken from
%   B. N. Cuffin and D. Cohen. Comparion of the Magnetoencephalogram and the Electroencephalogram. Electroencephalogr Clin Neurophysiol, 47:131-146, 1979.
%
% $Log: eeg_leadfield4_prepare.m,v $
% Revision 1.4  2003/07/29 16:02:07  roberto
% minor cosmetic change in the constants
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

r1 = vol.r(1); c1 = vol.c(1);
r2 = vol.r(2); c2 = vol.c(2);
r3 = vol.r(3); c3 = vol.c(3);
r4 = vol.r(4); c4 = vol.c(4);

if nargin==1
  Nmax = 60;
end

% these are the constants of cuffin1979
k1 = c1/c2;
k2 = c2/c3;
k3 = c3/c4;

for n=1:Nmax
  % according to lutkenhoner1992 the constant C is
  % lut_t(n) = ((n*c1/c2+n+1)*(n*c2/c3+n+1)+n*(n+1)*(c1/c2-1)*(c2/c3-1)*(r1/r2)^(2*n+1)) * ...
  %        ((n*c3/c4+n+1)+(n+1)*(c3/c4-1)*(r3/r4)^(2*n+1)) + ...
  %        ((c1/c2-1)*((n+1)*c2/c3+n)*(r1/r3)^(2*n+1)+(n*c1/c2+n+1)*(c2/c3-1)*(r2/r3)^(2*n+1)) * ...
  %        (n+1)*(n*(c3/c4-1)+((n+1)*c3/c4+n)*(r3/r4)^(2*n+1));
  % which can be rewritten as
  lut_t(n) = ((n*k1+n+1)*(n*k2+n+1)+n*(n+1)*(k1-1)*(k2-1)*(r1/r2)^(2*n+1)) * ...
         ((n*k3+n+1)+(n+1)*(k3-1)*(r3/r4)^(2*n+1)) + ...
         ((k1-1)*((n+1)*k2+n)*(r1/r3)^(2*n+1)+(n*k1+n+1)*(k2-1)*(r2/r3)^(2*n+1)) * ...
         (n+1)*(n*(k3-1)+((n+1)*k3+n)*(r3/r4)^(2*n+1));

end

% for debugging purposes, it can also give the constants of cuffin19979
if nargout>1
  % some extra constants of cuffin1979
  b = r1/r4;
  c = r2/r4;
  d = r3/r4;

  % according to cuffin1979 the constant Tau is (re-entered on 25 sept 2002)
  % but this requires also slightly other constants in the eeg_leadfield4 function
  for n=1:Nmax
    cuf_t(n) = d^(2*n+1) * (b^(2*n+1)*n*(k1-1)*(k2-1)*(n+1)...
           + c^(2*n+1)*(k1*n+n+1)*(k2*n+n+1))...
           *((k3*n+n+1)+(n+1)*(k3-1)*d^(2*n+1))...
           +(n+1)*c^(2*n+1)*(b^(2*n+1)*(k1-1)*(k2*n+k2+n)...
           +c^(2*n+1)*(k1*n+n+1)*(k2-1))...
           *(n*(k3-1)+(k3*n+k3+n)*d^(2*n+1));
  end
end
