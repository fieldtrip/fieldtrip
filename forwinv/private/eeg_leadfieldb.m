function [lf] = eeg_leadfieldb(pos, elc, vol)

% EEG_LEADFIELDB computes the electric leadfield for a dipole in a volume
% using the boundary element method
%
% [lf] = eeg_leadfieldb(pos, elc, vol)
%
% with the input arguments
%   pos		position dipole (1x3 or Nx3)
%   elc		position electrodes (optional, can be empty)
%   vol		volume conductor model
%
% the volume conductor model is a structure and should have the fields
%   vol.bnd	structure array with vertices and triangles of each boundary
%   vol.cond	conductivity of all compartments
%   vol.mat 	system matrix, which can include the electrode interpolation
%
% the compartment boundaries are described by a structure array with
%   vol.bnd(i).pnt
%   vol.bnd(i).pnt

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: eeg_leadfieldb.m,v $
% Revision 1.5  2009/04/23 15:06:14  roboos
% added patch from Cristiano
%
% Revision 1.4  2009/03/30 15:06:14  roboos
% added the patch from Alexandre to support openmeeg
%
% Revision 1.3  2009/02/02 13:12:53  roboos
% small fix, cpbem->bemcp
%
% Revision 1.2  2009/02/02 12:59:26  roboos
% added bemcp implementation, do not use vol.tra any more, made some changes to the checks on the input structures, use voltype in switch ladder
%
% Revision 1.1  2009/01/21 10:32:38  roboos
% moved from forwinv/* and forwinv/mex/* directory to forwinv/private/* to make the CVS layout consistent with the release version
%
% Revision 1.3  2008/04/14 20:54:35  roboos
% be more explicit about BEM type
%
% Revision 1.2  2005/12/06 11:41:07  roboos
% added support for dipoli models
% restructured the whole code
% combine electrode transfer and system matrix in a single matrix to speed forward computations up
%
% Revision 1.1  2003/10/07 08:40:22  roberto
% made separate function for BEM, based on part of eeg_leadfield.m
%


% do some sanity checks
if ~isfield(vol, 'bnd')
  error('there are no compartment boundaries present');
end

if length(vol.bnd)~=length(vol.cond)
  error('the number of compartments in the volume in ambiguous');
end

if ~isfield(vol, 'mat')
  error('there is no BEM system matrix present');
end

% determine the number of compartments
ncmp = length(vol.bnd);

% the number of rows in the leadfield matrix should either correspond to
% the number of electrodes, to the number of vertices of the skin
% compartment or to the total number of vertices
nelc  = size(elc, 1);
nskin = size(vol.bnd(vol.skin).pnt,1);
nall  = 0;
for i=1:ncmp
  nall = nall + size(vol.bnd(i).pnt,1);
end
if size(vol.mat,1)==nelc
  % the output leadfield corresponds to the number of electrodes
elseif size(vol.mat,1)==nskin
  % the output leadfield corresponds to the number skin vertices
elseif size(vol.mat,1)==nall
  % the output leadfield corresponds to the total number of vertices
elseif strcmp(voltype(vol),'openmeeg')
  % this is handled differently, although at the moment I don't know why
else
  error('unexpected size of vol.mat')
end

% determine the conductivity of the source compartment
cond = vol.cond(vol.source);

% compute the infinite medium potential on all vertices
switch voltype(vol)
  case 'avo'
    % the system matrix was computed using code from Adriaan van Oosterom
    % the code by Adriaan van Oosterom does not implement isolated source approach
    lf = [];
    for i=1:ncmp
      lf = [lf; inf_medium_leadfield(pos, vol.bnd(i).pnt, mean(vol.sigmas(i,:)))];
    end

  case 'dipoli'
    % the system matrix was computed using Thom Oostendorp's DIPOLI
    % concatenate the vertices of all compartment boundaries in a single Nx3 matrix
    pnt = [];
    for i=1:ncmp
      pnt = [pnt; vol.bnd(i).pnt];
    end
    % dipoli incorporates the conductivity into the system matrix
    lf = inf_medium_leadfield(pos, pnt, 1);

  case 'asa'
    % the system matrix was computed using ASA from www.ant-neuro.com
    % concatenate the vertices of all compartment boundaries in a single Nx3 matrix
    pnt = [];
    for i=1:ncmp
      pnt = [pnt; vol.bnd(i).pnt];
    end
    % assume that isolated potential approach was used
    lf = inf_medium_leadfield(pos, pnt, cond);

  case 'bemcp'
    % the system matrix was computed using code from Christopher Phillips
    cond = [vol.cond 0]; % add the conductivity of air for simplicity
    lf = cell(1,ncmp);
    % loop over boundaries and compute the leadfield for each
    for i=1:ncmp
      co = (cond(i)+cond(i+1))/2 ;
      lf{i} = inf_medium_leadfield(pos, vol.bnd(i).pnt, co);
    end
    % concatenate the leadfields
    lf = cat(1, lf{:});

   case 'openmeeg'
     % the system matrix is computed using OpenMEEG (Symmetric BEM)
     lf = openmeeg_lf_eeg(pos, elc, vol);

  otherwise
    error('unsupported type of volume conductor (%s)\n', voltype(vol));
end % switch voltype

if isfield(vol, 'mat') && ~voltype(vol, 'openmeeg')
  % compute the bounded medium potential on all vertices
  % this may include the bilinear interpolation from vertices towards electrodes
  lf = vol.mat * lf;
end

