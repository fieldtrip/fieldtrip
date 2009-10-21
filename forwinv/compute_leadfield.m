function [lf] = compute_leadfield(pos, sens, vol, varargin)

% COMPUTE_LEADFIELD computes a forward solution for a dipole in a a volume
% conductor model. The forward solution is expressed as the leadfield
% matrix (Nchan*3), where each column corresponds with the potential or field
% distributions on all sensors for one of the x,y,z-orientations of the
% dipole.
%
% Use as
%   [lf] = compute_leadfield(pos, sens, vol, ...)
% with input arguments
%   pos    position dipole (1x3 or Nx3)
%   sens   structure with gradiometer or electrode definition
%   vol    structure with volume conductor definition
%
% The vol structure represents a volume conductor model, its contents
% depend on the type of model. The sens structure represents a sensor
% arary, i.e. EEG electrodes or MEG gradiometers.
%
% It is possible to compute a simultaneous forward solution for EEG and MEG
% by specifying sens and grad as two cell-arrays, e.g.
%   sens = {senseeg, sensmeg}
%   vol  = {voleeg, volmeg }
% This results in the computation of the leadfield of the first element of
% sens and vol, followed by the second, etc. The leadfields of the
% different imaging modalities are concatenated.
%
% Additional input arguments can be specified as key-value pairs, supported
% optional arguments are
%   'reducerank'      = 'no' or number
%   'normalize'       = 'no', 'yes' or 'column'
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%
% Depending on the specific input arguments for the sensor and volume, this
% function will select the appropriate low-level EEG or MEG forward model.
% The leadfield matrix for EEG will have an average reference over all the
% electrodes.
%
% The supported forward solutions for MEG are
%   single sphere (Cuffin and Cohen, 1977)
%   multiple spheres with one sphere per channel (Huang et al, 1999)
%   realistic single shell using superposition of basis functions (Nolte, 2003)
%   using leadfield interpolation on precomputed grid
%   Boundary Element Method (BEM) using Neuromag meg_calc toolbox
%
% The supported forward solutions for EEG are
%   single sphere
%   multiple concentric spheres (max. 4)
%   using leadfield interpolation on precomputed grid
%   Boundary Element Method (BEM) using ASA to precompute the sytem matrix
%   Boundary Element Method (BEM) using Neuromag meg_calc toolbox
%
% References to implemented methods:
%   Cuffin BN, Cohen D.
%   Magnetic fields of a dipole in special volume conductor shapes
%   IEEE Trans Biomed Eng. 1977 Jul;24(4):372-81.
%
%   Nolte G.
%   The magnetic lead field theorem in the quasi-static approximation and its use for magnetoencephalography forward calculation in realistic volume conductors
%   Phys Med Biol. 2003 Nov 21;48(22):3637-52
%
%   Huang MX, Mosher JC, Leahy RM.
%   A sensor-weighted overlapping-sphere head model and exhaustive head model comparison for MEG
%   Phys Med Biol. 1999 Feb;44(2):423-40

% Copyright (C) 2004-2008, Robert Oostenveld
%
% $Log: compute_leadfield.m,v $
% Revision 1.30  2009/05/25 08:17:15  roboos
% small change in consistency test for multisphere MEG volume conductor
%
% Revision 1.29  2009/03/30 15:06:14  roboos
% added the patch from Alexandre to support openmeeg
%
% Revision 1.28  2009/02/02 13:05:44  roboos
% addec bemcp
% give warning once in case of eeg infinite medium
%
% Revision 1.27  2008/07/22 10:17:15  roboos
% replaced identical with strcmp
%
% Revision 1.26  2008/07/21 20:28:44  roboos
% added check on units (mm/cm/m) of the sensor array and volume conductor, give error if inconsistent
%
% Revision 1.25  2008/05/13 19:24:11  roboos
% consistently removed support for pnt1/pnt2 gradiometer description
%
% Revision 1.24  2008/04/30 13:47:20  roboos
% removed support for pnt1+pnt2
%
% Revision 1.23  2008/04/15 20:36:21  roboos
% added explicit handling of various BEM implementations, i.e. for all voltype variants
%
% Revision 1.22  2008/04/11 13:16:16  roboos
% added support for simultaneous EEG and MEG
%
% Revision 1.21  2008/03/18 13:20:42  roboos
% updated documentation and some error messages
%
% Revision 1.20  2008/03/05 16:24:30  roboos
% use the voltype helper function
% added option for 'normalizeparam'
%
% Revision 1.19  2007/07/25 08:32:25  roboos
% switched to using senstype helper function
%
% Revision 1.18  2006/10/16 15:21:04  roboos
% small change in a comment
%
% Revision 1.17  2006/10/12 08:57:34  roboos
% added support for multiple dipoles for which the positions are in one long 1x(3*N) vector
% changed the functino call to neuromag megfield (not verified)
% renamed some internal variables
% only determine the number of dipoles (Ndipoles) once, and not in every section again
%
% Revision 1.16  2006/10/04 08:12:40  roboos
% changed some '&' into '&&', only apply reducerank when smaller than number of leadfield components (for efficiency)
%
% Revision 1.15  2006/09/07 12:43:22  roboos
% added code for multisphere, contributed by Punita Christopher (based on BrainStorm)
%
% Revision 1.14  2006/03/21 11:24:16  roboos
% multiply the nolte leadfield with sens.tra to transform from coils to gradiometers
%
% Revision 1.13  2006/03/21 09:40:08  roboos
% implemented support for Guido Nolte's method directly using his code
% improved documentation, added literature references
%
% Revision 1.12  2006/02/28 13:35:14  roboos
% added column-wise normalization for leadfield
%
% Revision 1.11  2006/02/14 09:42:31  roboos
% added a snippet of code to support forward computations using the Neuromag meg-calc toolbox
%
% Revision 1.10  2006/01/20 09:41:49  roboos
% changed assignment of sens.pnt (add structure field to existing double) according to matlab 7.1. recommendation (prevents warning)
%
% Revision 1.9  2005/11/16 09:56:17  roboos
% added a semicolon to prevent output on screen
%
% Revision 1.8  2005/11/08 11:05:14  roboos
% added support for normalize and reducerank as additional options (using key-value varargin)
%
% Revision 1.7  2005/07/29 07:19:47  roboos
% improved detection of meg data (ismeg if has ori and pos)
%
% Revision 1.6  2005/06/08 16:35:15  roboos
% corrected the check for the number of spheres in the grad and volume
%
% Revision 1.5  2005/02/21 08:01:19  roboos
% removed spaces at the end of line, no code changes
%
% Revision 1.4  2005/02/08 11:59:14  roboos
% added an extra check on the input
%
% Revision 1.3  2004/10/25 16:20:13  roboos
% fixed bug in selection of electrode positions from sens structure
%
% Revision 1.2  2004/09/21 15:31:51  roboos
% renamed grad into sens
%
% Revision 1.1  2004/08/19 08:18:20  roboos
% new implementations, generalized for EEG and MEG
%

persistent warning_issued;

if iscell(sens) && iscell(vol) && numel(sens)==numel(vol)
  % this represents combined EEG and MEG sensors, where each modality has its own volume conduction model
  lf = cell(1,numel(sens));
  for i=1:length(sens)
    lf{i} = compute_leadfield(pos, sens{i}, vol{i}, varargin{:});
  end
  lf = cat(1, lf{:});
  return;
end

if ~isstruct(sens) && size(sens,2)==3
  % definition of electrode positions only, restructure it
  sens = struct('pnt', sens);
end

% determine whether it is EEG or MEG
iseeg = senstype(sens, 'eeg');
ismeg = senstype(sens, 'meg');

% get the optional input arguments
reducerank     = keyval('reducerank', varargin); if isempty(reducerank), reducerank = 'no'; end
normalize      = keyval('normalize' , varargin); if isempty(normalize ), normalize  = 'no'; end
normalizeparam = keyval('normalizeparam', varargin); if isempty(normalizeparam ), normalizeparam = 0.5; end

% multiple dipoles can be represented either as a 1x(N*3) vector or as a
% as a Nx3 matrix, i.e. [x1 y1 z1 x2 y2 z2] or [x1 y1 z1; x2 y2 z2]
Ndipoles = numel(pos)/3;
if all(size(pos)==[1 3*Ndipoles])
  pos = reshape(pos, 3, Ndipoles)';
end

if isfield(vol, 'unit') && isfield(sens, 'unit') && ~strcmp(vol.unit, sens.unit)
  error('inconsistency in the units of the volume conductor and the sensor array');
end

if ismeg && iseeg
  % this is something that could be implemented relatively easily
  error('simultaneous EEG and MEG not supported');

elseif ~ismeg && ~iseeg
  error('the input does not look like EEG, nor like MEG');

elseif ismeg
  switch voltype(vol)

    case 'singlesphere'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % MEG single-sphere volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      pnt = sens.pnt; % position of each coil
      ori = sens.ori; % orientation of each coil

      if isfield(vol, 'o')
        % shift dipole and magnetometers to origin of sphere
        pos = pos - repmat(vol.o, Ndipoles, 1);
        pnt = pnt - repmat(vol.o, size(pnt,1), 1);
      end

      if Ndipoles>1
        % loop over multiple dipoles
        lf = zeros(size(pnt,1),3*Ndipoles);
        for i=1:Ndipoles
          lf(:,(3*i-2):(3*i)) = meg_leadfield1(pos(i,:), pnt, ori);
        end
      else
        % only single dipole
        lf = meg_leadfield1(pos, pnt, ori);
      end

      if isfield(sens, 'tra')
        % this appears to be the modern complex gradiometer definition
        % construct the channels from a linear combination of all magnetometers
        lf = sens.tra * lf;
      end

    case 'multisphere'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % MEG multi-sphere volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ncoils = length(sens.pnt);

      if size(vol.r, 1)~=ncoils
        error('number of spheres is not equal to the number of coils')
      end

      if size(vol.o, 1)~=ncoils
        error('number of spheres is not equal to the number of coils');
      end

      lf = zeros(ncoils, 3*Ndipoles);
      for chan=1:ncoils
        for dip=1:Ndipoles
          % shift dipole and magnetometer coil to origin of sphere
          dippos = pos(dip,:)       - vol.o(chan,:);
          chnpos = sens.pnt(chan,:) - vol.o(chan,:);
          tmp = meg_leadfield1(dippos, chnpos, sens.ori(chan,:));
          lf(chan,(3*dip-2):(3*dip)) = tmp;
        end
      end

      if isfield(sens, 'tra')
        % this appears to be the modern complex gradiometer definition
        % construct the channels from a linear combination of all magnetometers
        lf = sens.tra * lf;
      end

    case 'neuromag'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % use external Neuromag toolbox for forward computation
      % this requires that "megmodel" is initialized, which is done in PREPARE_VOL_SENS
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % compute the forward model for all channels
      % tmp1 = ones(1, Ndipoles);
      % tmp2 = 0.01*pos';  %convert to cm
      % lf = megfield([tmp2 tmp2 tmp2],[[1 0 0]'*tmp1 [0 1 0]'*tmp1 [0 0 1]'*tmp1]);
      for dip=1:Ndipoles
        R = 0.01*pos(i,:)'; % convert from cm to m
        Qx = [1 0 0];
        Qy = [0 1 0];
        Qz = [0 0 1];
        lf(:,(3*(dip-1)+1)) = megfield(R, Qx);
        lf(:,(3*(dip-1)+2)) = megfield(R, Qy);
        lf(:,(3*(dip-1)+3)) = megfield(R, Qz);
      end
      % select only those channels from the forward model that are part of the gradiometer definition
      lf = lf(vol.chansel,:);

    case 'nolte'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % use code from Guido Nolte for the forward computation
      % this requires that "meg_ini" is initialized, which is done in PREPARE_VOL_SENS
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % the dipole position and orientation should be combined in a single matrix
      % furthermore, here I want to compute the leadfield for each of the
      % orthogonzl x/y/z directions
      dippar = zeros(Ndipoles*3, 6);
      for i=1:Ndipoles
        dippar((i-1)*3+1,:) = [pos(i,:) 1 0 0];  % single dipole, x-orientation
        dippar((i-1)*3+2,:) = [pos(i,:) 0 1 0];  % single dipole, y-orientation
        dippar((i-1)*3+3,:) = [pos(i,:) 0 0 1];  % single dipole, z-orientation
      end
      % compute the leadfield for each individual coil
      lf = meg_forward(dippar,vol.forwpar);
      if isfield(sens, 'tra')
        % compute the leadfield for each gradiometer (linear combination of coils)
        lf = sens.tra * lf;
      end

    case 'infinite'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % magnetic dipole instead of electric (current) dipole in an infinite vacuum
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if isempty(warning_issued)
        % give the warning only once
        warning('assuming magnetic dipole in an infinite vacuum');
        warning_issued = 1;
      end

      pnt = sens.pnt; % position of each coil
      ori = sens.ori; % orientation of each coil

      if Ndipoles>1
        % loop over multiple dipoles
        lf = zeros(size(pnt,1),3*Ndipoles);
        for i=1:Ndipoles
          lf(:,(3*i-2):(3*i)) = magnetic_dipole(pos(i,:), pnt, ori);
        end
      else
        % only single dipole
        lf = magnetic_dipole(pos, pnt, ori);
      end

      if isfield(sens, 'tra')
        % construct the channels from a linear combination of all magnetometer coils
        lf = sens.tra * lf;
      end

    otherwise
      error('unsupported volume conductor model for MEG');
  end % switch voltype for MEG

elseif iseeg
  switch voltype(vol)

    case 'multisphere'
      % Based on the approximation of the potential due to a single dipole in
      % a multishell sphere by three dipoles in a homogeneous sphere, code
      % contributed by Punita Christopher

      Nelec = size(sens.pnt,1);
      Nspheres = length(vol.r);

      % the center of the spherical volume conduction model does not have
      % to be in the origin, therefore shift the spheres, the electrodes
      % and the dipole
      if isfield(vol, 'o')
        center = vol.o;
      else
        center = [0 0 0];
      end

      % sort the spheres from the smallest to the largest
      % furthermore, the radius should be one (?)
      [radii, indx] = sort(vol.r/max(vol.r));
      sigma = vol.c(indx);
      r   = (sens.pnt-repmat(center, Nelec, 1))./max(vol.r);
      pos = pos./max(vol.r);

      if Ndipoles>1
        % loop over multiple dipoles
        lf = zeros(Nelec,3*Ndipoles);
        for i=1:Ndipoles
          rq = pos(i,:) - center;
          % compute the potential for each dipole ortientation
          % it would be much more efficient to change the punita function
          q1 = [1 0 0]; lf(:,(3*i-2)) = multisphere(Nspheres, radii, sigma, r, rq, q1);
          q1 = [0 1 0]; lf(:,(3*i-1)) = multisphere(Nspheres, radii, sigma, r, rq, q1);
          q1 = [0 0 1]; lf(:,(3*i  )) = multisphere(Nspheres, radii, sigma, r, rq, q1);
        end
      else
        % only single dipole
        lf = zeros(Nelec,3);
        rq = pos - center;
        % compute the potential for each dipole ortientation
        % it would be much more efficient to change the punita function
        q1 = [1 0 0] ; lf(:,1) = multisphere(Nspheres, radii, sigma, r, rq, q1);
        q1 = [0 1 0] ; lf(:,2) = multisphere(Nspheres, radii, sigma, r, rq, q1);
        q1 = [0 0 1] ; lf(:,3) = multisphere(Nspheres, radii, sigma, r, rq, q1);
      end

    case {'singlesphere', 'concentric'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % EEG spherical volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % FIXME, this is not consistent between spherical and BEM
      % sort the spheres from the smallest to the largest
      [vol.r, indx] = sort(vol.r);
      vol.c = vol.c(indx);

      Nspheres = length(vol.c);
      if length(vol.r)~=Nspheres
        error('the number of spheres in the volume conductor model is ambiguous');
      end

      if isfield(vol, 'o')
        % shift the origin of the spheres, electrodes and dipole
        sens.pnt = sens.pnt - repmat(vol.o, size(sens.pnt,1), 1);
        pos = pos - repmat(vol.o, Ndipoles, 1);
      end

      if Nspheres==1
        if Ndipoles>1
          % loop over multiple dipoles
          lf = zeros(size(sens.pnt,1),3*Ndipoles);
          for i=1:Ndipoles
            lf(:,(3*i-2):(3*i)) = eeg_leadfield1(pos(i,:), sens.pnt, vol);
          end
        else
          % only single dipole
          lf = eeg_leadfield1(pos, sens.pnt, vol);
        end

      elseif Nspheres==2
        vol.r = [vol.r(1) vol.r(2) vol.r(2) vol.r(2)];
        vol.c = [vol.c(1) vol.c(2) vol.c(2) vol.c(2)];
        if Ndipoles>1
          % loop over multiple dipoles
          lf = zeros(size(sens.pnt,1),3*Ndipoles);
          for i=1:Ndipoles
            lf(:,(3*i-2):(3*i)) = eeg_leadfield4(pos(i,:), sens.pnt, vol);
          end
        else
          % only single dipole
          lf = eeg_leadfield4(pos, sens.pnt, vol);
        end

      elseif Nspheres==3
        vol.r = [vol.r(1) vol.r(2) vol.r(3) vol.r(3)];
        vol.c = [vol.c(1) vol.c(2) vol.c(3) vol.c(3)];
        if Ndipoles>1
          % loop over multiple dipoles
          lf = zeros(size(sens.pnt,1),3*Ndipoles);
          for i=1:Ndipoles
            lf(:,(3*i-2):(3*i)) = eeg_leadfield4(pos(i,:), sens.pnt, vol);
          end
        else
          % only single dipole
          lf = eeg_leadfield4(pos, sens.pnt, vol);
        end

      elseif Nspheres==4
        if Ndipoles>1
          % loop over multiple dipoles
          lf = zeros(size(sens.pnt,1),3*Ndipoles);
          for i=1:Ndipoles
            lf(:,(3*i-2):(3*i)) = eeg_leadfield4(pos(i,:), sens.pnt, vol);
          end
        else
          % only single dipole
          lf = eeg_leadfield4(pos, sens.pnt, vol);
        end

      else
        error('more than 4 concentric spheres are not supported');
      end

    case {'bem', 'dipoli', 'asa', 'avo', 'bemcp', 'openmeeg'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % EEG boundary element method volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      lf = eeg_leadfieldb(pos, sens.pnt, vol);

    case 'infinite'
      % the conductivity of the medium is not known
      if isempty(warning_issued)
        % give the warning only once
        warning('assuming electric dipole in an infinite medium with unit conductivity');
        warning_issued = 1;
      end
      lf = inf_medium_leadfield(pos, sens.pnt, 1);

    otherwise
      error('unsupported volume conductor model for EEG');
  end % switch voltype for EEG

  % compute average reference for EEG leadfield
  avg = mean(lf, 1);
  lf  = lf - repmat(avg, size(lf,1), 1);

end % iseeg or ismeg

% optionally apply leadfield rank reduction
if ~strcmp(reducerank, 'no') && reducerank<size(lf,2)
  % decompose the leadfield
  [u, s, v] = svd(lf);
  r = diag(s);
  s(:) = 0;
  for j=1:reducerank
    s(j,j) = r(j);
  end
  % recompose the leadfield with reduced rank
  lf = u * s * v';
end

% optionally apply leadfield normaliziation
if strcmp(normalize, 'yes')
  if normalizeparam==0.5
    % normalize the leadfield by the Frobenius norm of the matrix
    % this is the same as below in case normalizeparam is 0.5
    nrm = norm(lf, 'fro');
  else
    % normalize the leadfield by sum of squares of the elements of the leadfield matrix to the power "normalizeparam"
    % this is the same as the Frobenius norm if normalizeparam is 0.5
    nrm = sum(lf(:).^2)^normalizeparam;
  end
  if nrm>0
    lf = lf ./ nrm;
  end
elseif strcmp(normalize, 'column')
  % normalize each column of the leadfield by its norm
  for j=1:size(lf,2)
    nrm = sum(lf(:,j).^2)^normalizeparam;
    lf(:,j) = lf(:,j)./nrm;
  end
end
