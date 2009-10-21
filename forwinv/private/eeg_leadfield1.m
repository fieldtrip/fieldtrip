function [lf, lforig] = eeg_leadfield1(R, elc, vol);

% EEG_LEADFIELD1 electric leadfield for a dipole in a single sphere
%
% [lf] = eeg_leadfield1(R, elc, vol)
%
% with input arguments
%   R		position dipole (vector of length 3)
%   elc		position electrodes
% and vol being a structure with the elements
%   vol.r	radius of sphere
%   vol.c	conductivity of sphere

% Copyright (C) 2002, Robert Oostenveld
%
% this implementation is adapted from
%   Luetkenhoener, Habilschrift '92
% the original reference is
%   R. Kavanagh, T. M. Darccey, D. Lehmann, and D. H. Fender. Evaluation of methods for three-dimensional localization of electric sources in the human brain. IEEE Trans Biomed Eng, 25:421-429, 1978.
%
% $Log: eeg_leadfield1.m,v $
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.11  2008/12/24 13:33:45  roboos
% fixed typo in documentation
%
% Revision 1.10  2006/02/09 08:30:28  roboos
% scale with conductivity in case of dipole in the origin (thanks to Denise van Barneveld)
%
% Revision 1.9  2006/01/20 09:40:09  roboos
% disabled the dipole-inside-brain check (for EEGLAB)
% use precomputed c5 to capture divide by zero situation
%
% Revision 1.8  2004/10/25 16:25:25  roboos
% re-enabled check for dipole outside head
%
% Revision 1.7  2003/12/05 09:48:45  roberto
% disabled error check for dipole outside brain [for EEGLAB]
%
% Revision 1.6  2003/12/04 10:43:51  roberto
% added error when dipole is outside brain compartment
%
% Revision 1.5  2003/07/29 16:04:55  roberto
% fixed bug in determining whether electrodes are lying on sphere surface
%
% Revision 1.4  2003/07/29 15:53:09  roberto
% minor change to the projection of the electrodes to the sphere (now default)
%
% Revision 1.3  2003/07/29 15:42:03  roberto
% found and fixed a bug in the handling of electrodes and dipoles on one line
% changed the name of the constants and added a constant
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

Nchans = size(elc, 1);
lf = zeros(Nchans,3);

% always take the outermost sphere, this makes comparison with the 4-sphere computation easier
[vol.r, indx] = max(vol.r);
vol.c = vol.c(indx);

% check whether the electrode ly on the sphere, allowing 0.5% tolerance
dist = sqrt(sum(elc.^2,2));
if any(abs(dist-vol.r)>vol.r*0.005)
  warning('electrodes do not ly on sphere surface -> using projection')
end
elc = vol.r * elc ./ [dist dist dist];

% check whether the dipole is inside the brain [disabled for EEGLAB]
% if sqrt(sum(R.^2))>=vol.r
%   error('dipole is outside the brain compartment');
% end

c0 = norm(R);
c1 = vol.r;
c2 = 4*pi*c0^2*vol.c;

if c0==0
  % the dipole is in the origin, this can and should be handeled as an exception
  [phi, el] = cart2sph(elc(:,1), elc(:,2), elc(:,3));
  theta = pi/2 - el;
  lf(:,1) = sin(theta).*cos(phi);
  lf(:,2) = sin(theta).*sin(phi);
  lf(:,3) = cos(theta);
  % the potential in a homogenous sphere is three times the infinite medium potential
  lf = 3/(c1^2*4*pi*vol.c)*lf;

else
  for i=1:Nchans
    % use another name for the electrode, in accordance with lutkenhoner1992
    r = elc(i,:);

    c3 = r-R;
    c4 = norm(c3);
    c5 = c1^2 * c0^2 - dot(r,R)^2;   % lutkenhoner A.11
    c6 = c0^2*r - dot(r,R)*R;        % lutkenhoner, just after A.17

    % the original code reads (cf. lutkenhoner1992 equation A.17)
    % lf(i,:) = ((dot(R, r/norm(r) - (r-R)/norm(r-R))/(norm(cross(r,R))^2) + 2/(norm(r-R)^3)) * cross(R, cross(r, R)) + ((norm(r)^2-norm(R)^2)/(norm(r-R)^3) - 1/norm(r)) * R) / (4*pi*vol.c(1)*norm(R)^2);

    % but more efficient execution of the code is achieved by some precomputations
    if c5<1000*eps
      % the dipole lies on a single line with the electrode
      lf(i,:) = (2/c4^3 * c6 + ((c1^2-c0^2)/c4^3 - 1/c1) * R) / c2;
    else
      % nothing wrong, do the complete computation
      lf(i,:) = ((dot(R, r/c1 - c3/c4)/c5 + 2/c4^3) * c6 + ((c1^2-c0^2)/c4^3 - 1/c1) * R) / c2;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product
function [c] = cross(a,b)
c = [a(2)*b(3)-a(3)*b(2) a(3)*b(1)-a(1)*b(3) a(1)*b(2)-a(2)*b(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast dot product
function [c] = dot(a,b)
c = sum(a.*b);

