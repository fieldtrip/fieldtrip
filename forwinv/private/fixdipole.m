function dip = fixdipole(dip)

% FIXDIPOLE ensures that the dipole position and moment are
% consistently represented throughout FieldTrip functions.

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: fixdipole.m,v $
% Revision 1.1  2009/07/02 15:35:55  roboos
% helper function for consistent dipole representation
%

[m, n] = size(dip.pos);

if n==3
  % the input representation is Nx3, which is what we want
elseif m==3
  % it is possible to translate it into a Nx3 unambiguously
  warning('input dipole positions should be specified as Nx3 matrix');
  dip.pos = dip.pos';
elseif m==1
  % it is possible to translate it into a Nx3 unambiguously
  warning('input dipole positions should be specified as Nx3 matrix');
  dip.pos = reshape(dip.pos, 3, n/3)';
else
  % it is not clear how to convert to a Nx3 matrix
  error('input dipole positions should be specified as Nx3 matrix');
end

if isfield(dip, 'mom')
  ndip = size(dip.pos,1);
  if numel(dip.mom)==ndip*3
    ntime = 1;
  else
    ntime = numel(dip.mom)/(ndip*3);
  end
  dip.mom = reshape(dip.mom, ndip*3, ntime);
end

