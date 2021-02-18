function [lf] = ft_compute_leadfield(dippos, sens, headmodel, varargin)

% FT_COMPUTE_LEADFIELD computes a forward solution for a dipole in a a volume
% conductor model. The forward solution is expressed as the leadfield matrix
% (Nchan*3), where each column corresponds with the potential or field distributions
% on all sensors for one of the x,y,z-orientations of the dipole.
%
% Use as
%   [lf] = ft_compute_leadfield(dippos, sens, headmodel, ...)
% with input arguments
%   dippos    = position dipole (1*3 or Ndip*3)
%   sens      = structure with gradiometer or electrode definition
%   headmodel = structure with volume conductor definition
%
% The headmodel represents a volume conductor model, its contents
% depend on the type of model. The sens structure represents a sensor
% array, i.e. EEG electrodes or MEG gradiometers.
%
% It is possible to compute a simultaneous forward solution for EEG and MEG
% by specifying sens and grad as two cell-arrays, e.g.
%   sens       = {senseeg, sensmeg}
%   headmodel  = {voleeg,  volmeg}
% This results in the computation of the leadfield of the first element of
% sens and headmodel, followed by the second, etc. The leadfields of the
% different imaging modalities are subsequently concatenated.
%
% Additional input arguments can be specified as key-value pairs, supported
% optional arguments are
%   'reducerank'      = 'no' or number (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% The leadfield weight may be used to specify a (normalized) corresponding surface
% area for each dipole, e.g. when the dipoles represent a folded cortical surface
% with varying triangle size.
%
% Depending on the specific input arguments for the sensor and volume, this function
% will select the appropriate low-level EEG or MEG forward model. The leadfield
% matrix for EEG will have an average reference over all the electrodes.
%
% The supported forward solutions for MEG are
%   infinite homogenous medium
%   single sphere (Cuffin and Cohen, 1977)
%   multiple spheres with one sphere per channel (Huang et al, 1999)
%   realistic single shell using superposition of basis functions (Nolte, 2003)
%   leadfield interpolation using a precomputed sourcemodel
%   boundary element method (BEM)
%
% The supported forward solutions for EEG are
%   infinite homogenous medium
%   infinite halfspace homogenous medium
%   single sphere
%   multiple concentric spheres (up to 4 spheres)
%   leadfield interpolation using a precomputed sourcemodel
%   boundary element method (BEM)
%   finite element method (FEM)
%
% See also FT_PREPARE_VOL_SENS, FT_HEADMODEL_ASA, FT_HEADMODEL_BEMCP,
% FT_HEADMODEL_CONCENTRICSPHERES, FT_HEADMODEL_DIPOLI, FT_HEADMODEL_HALFSPACE,
% FT_HEADMODEL_INFINITE, FT_HEADMODEL_LOCALSPHERES, FT_HEADMODEL_OPENMEEG,
% FT_HEADMODEL_SINGLESHELL, FT_HEADMODEL_SINGLESPHERE,
% FT_HEADMODEL_HALFSPACE, FT_HEADMODEL_DUNEURO

% Copyright (C) 2004-2020, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if iscell(sens) && iscell(headmodel) && numel(sens)==numel(headmodel)
  % this represents combined EEG, ECoG and/or MEG
  % use recursion to compute all leadfields
  lf = cell(1, numel(sens));
  for i=1:length(sens)
    lf{i} = ft_compute_leadfield(dippos, sens{i}, headmodel{i}, varargin{:});
  end
  lf = cat(1, lf{:});
  return;
end

% get the optional input arguments
reducerank      = ft_getopt(varargin, 'reducerank'); % default is handled below
backproject     = ft_getopt(varargin, 'backproject', 'yes');
normalize       = ft_getopt(varargin, 'normalize' , 'no');
normalizeparam  = ft_getopt(varargin, 'normalizeparam', 0.5);
weight          = ft_getopt(varargin, 'weight');
chanunit        = ft_getopt(varargin, 'chanunit');   % this is something like V, T, or T/m
dipoleunit      = ft_getopt(varargin, 'dipoleunit'); % this is something like nA*m

if any(strcmp(varargin(1:2:end), 'unit'))
  ft_error('the ''unit'' option is not supported any more, please use ''chanunit''');
end
if any(strcmp(varargin(1:2:end), 'units'))
  ft_error('the ''units'' option is not supported any more, please use ''chanunit''');
end

if ~isstruct(sens) && size(sens, 2)==3
  % definition of electrode positions only, restructure it
  sens = struct('elecpos', sens);
end

% ft_prepare_vol_sens should be called prior to ft_compute_leadfield
% to ensure that the sens and headmodel are up to date, since the backward
% compatibility check should not be performed for each dipole location
% sens       = ft_datatype_sens(sens);
% headmodel  = ft_datatype_headmodel(headmodel);

% determine whether it is EEG or MEG
iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');

% determine the default for this option
if isempty(reducerank)
  if iseeg
    reducerank = 'no';    % for EEG
  elseif ismeg && ft_headmodeltype(headmodel, 'infinite')
    reducerank = 'no';    % for MEG with a magnetic dipole, e.g. a HPI coil
  else
    reducerank = 'yes';   % for MEG with a current dipole in a volume conductor
  end
end

% multiple dipoles can be represented either as a 1x(N*3) vector or as a
% as a Nx3 matrix, i.e. [x1 y1 z1 x2 y2 z2] or [x1 y1 z1; x2 y2 z2]
Ndipoles = numel(dippos)/3;
if all(size(dippos)==[1 3*Ndipoles])
  dippos = reshape(dippos, 3, Ndipoles)';
end

if isscalar(weight)
  weight = weight * ones(Ndipoles, 1);
end

if isfield(headmodel, 'unit') && isfield(sens, 'unit') && ~strcmp(headmodel.unit, sens.unit)
  ft_error('inconsistency in the units of the volume conductor and the sensor array');
end

if ismeg && iseeg
  % this is something that could be implemented relatively easily
  ft_error('simultaneous EEG and MEG not supported');

elseif ~ismeg && ~iseeg
  ft_error('the input does not look like EEG, nor like MEG');

elseif ismeg
  switch ft_headmodeltype(headmodel)

    case 'singlesphere'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % MEG single-sphere volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      coilpos = sens.coilpos; % position of each coil
      coilori = sens.coilori; % orientation of each coil

      if isfield(headmodel, 'o')
        % shift dipole and magnetometers to origin of sphere
        dippos  = dippos  - repmat(headmodel.o, Ndipoles, 1);
        coilpos = coilpos - repmat(headmodel.o, size(coilpos, 1), 1);
      end

      if Ndipoles>1
        % loop over multiple dipoles
        lf = zeros(size(coilpos, 1), 3*Ndipoles);
        for i=1:Ndipoles
          lf(:, (3*i-2):(3*i)) = meg_leadfield1(dippos(i, :), coilpos, coilori);
        end
      else
        % only single dipole
        lf = meg_leadfield1(dippos, coilpos, coilori);
      end

      if isfield(sens, 'tra')
        % this appears to be the modern complex gradiometer definition
        % construct the channels from a linear combination of all magnetometers
        lf = sens.tra * lf;
      end

    case 'localspheres'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % MEG multiple overlapping sphere volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ncoils = length(sens.coilpos);

      if size(headmodel.r, 1)~=ncoils
        ft_error('number of spheres is not equal to the number of coils')
      end

      if size(headmodel.o, 1)~=ncoils
        ft_error('number of spheres is not equal to the number of coils');
      end

      lf = zeros(ncoils, 3*Ndipoles);
      for coil=1:ncoils
        for dip=1:Ndipoles
          % shift dipole and magnetometer coil to origin of sphere
          tmppos  = dippos(dip, :) - headmodel.o(coil, :);
          coilpos = sens.coilpos(coil, :) - headmodel.o(coil, :);
          tmp = meg_leadfield1(tmppos, coilpos, sens.coilori(coil, :));
          lf(coil, (3*dip-2):(3*dip)) = tmp;
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
      % tmp2 = 0.01*dippos'; %convert to cm
      % lf = megfield([tmp2 tmp2 tmp2], [[1 0 0]'*tmp1 [0 1 0]'*tmp1 [0 0 1]'*tmp1]);
      for dip=1:Ndipoles
        R = 0.01*dippos(dip, :)'; % convert from cm to m
        Qx = [1 0 0];
        Qy = [0 1 0];
        Qz = [0 0 1];
        lf(:, (3*(dip-1)+1)) = megfield(R, Qx);
        lf(:, (3*(dip-1)+2)) = megfield(R, Qy);
        lf(:, (3*(dip-1)+3)) = megfield(R, Qz);
      end
      % select only those channels from the forward model that are part of the gradiometer definition
      lf = lf(headmodel.chansel, :);

    case 'singleshell'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % use code from Guido Nolte for the forward computation
      % this requires that "meg_ini" is initialized, which is done in PREPARE_VOL_SENS
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % the dipole position and orientation should be combined in a single matrix
      % furthermore, here I want to compute the leadfield for each of the
      % orthogonal x/y/z directions
      dippar = zeros(Ndipoles*3, 6);
      for i=1:Ndipoles
        dippar((i-1)*3+1, :) = [headmodel.forwpar.scale*dippos(i, :) 1 0 0]; % single dipole with unit strength, x-orientation
        dippar((i-1)*3+2, :) = [headmodel.forwpar.scale*dippos(i, :) 0 1 0]; % single dipole with unit strength, y-orientation
        dippar((i-1)*3+3, :) = [headmodel.forwpar.scale*dippos(i, :) 0 0 1]; % single dipole with unit strength, z-orientation
      end
      % compute the leadfield for each individual coil
      lf = meg_forward(dippar, headmodel.forwpar);
      % the leadfield is computed for cm units, convert it to the desired units
      lf = lf*headmodel.forwpar.scale^2;
      if isfield(sens, 'tra')
        % compute the leadfield for each gradiometer (linear combination of coils)
        lf = sens.tra * lf;
      end

    case 'openmeeg'
      ft_hastoolbox('openmeeg', 1);

      dsm         = ft_getopt(varargin, 'dsm');
      nonadaptive = ft_getopt(varargin, 'nonadaptive');

      [h2sens,ds2sens] = ft_sensinterp_openmeeg(dippos, headmodel, sens);
      if isempty(dsm)
        dsm            = ft_sysmat_openmeeg(dippos, headmodel, sens, nonadaptive);
      end
      lf               = ds2sens + h2sens*headmodel.mat*dsm;
      
      if isfield(sens, 'tra')
        % compute the leadfield for each gradiometer (linear combination of coils)
        lf = sens.tra * lf;
      end

    case {'infinite_magneticdipole', 'infinite'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % magnetic dipole instead of electric (current) dipole in an infinite vacuum
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      coilpos = sens.coilpos; % position of each coil
      coilori = sens.coilori; % orientation of each coil

      if Ndipoles>1
        % loop over multiple dipoles
        lf = zeros(size(coilpos, 1), 3*Ndipoles);
        for i=1:Ndipoles
          lf(:, (3*i-2):(3*i)) = magnetic_dipole(dippos(i, :), coilpos, coilori);
        end
      else
        % only single dipole
        lf = magnetic_dipole(dippos, coilpos, coilori);
      end

      if isfield(sens, 'tra')
        % construct the channels from a linear combination of all magnetometer coils
        lf = sens.tra * lf;
      end

    case {'infinite_currentdipole'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % current dipole in an infinite homogenous conducting medium
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      coilpos = sens.coilpos; % position of each coil
      coilori = sens.coilori; % orientation of each coil

      if Ndipoles>1
        % loop over multiple dipoles
        lf = zeros(size(coilpos, 1), 3*Ndipoles);
        for i=1:Ndipoles
          lf(:, (3*i-2):(3*i)) = current_dipole(dippos(i, :), coilpos, coilori);
        end
      else
        % only single dipole
        lf = current_dipole(dippos, coilpos, coilori);
      end

      if isfield(sens, 'tra')
        % construct the channels from a linear combination of all magnetometer coils
        lf = sens.tra * lf;
      end

    case {'duneuro'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % finite element method as implemented in software duneuro
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
      %TODO: involve unit checking -> not at this level

      % compute secondary leadfield numerically
      lf = leadfield_duneuro(dippos, headmodel, 'meg');

      % compute primary B-field analytically
      % permeability constant mu in si units
      mu = 4*pi*1e-7; %unit: Tm/A
      
      index = repmat(1:size(dippos,1),3,1);
      index = index(:);
      dipoles = [dippos(index,:) repmat(eye(3),size(dippos,1),1)];
      Bp = compute_B_primary(sens.coilpos, dipoles, sens.coilori);

      % compute full B-field
      lf = mu/(4*pi) * (Bp - lf);

      if isfield(sens, 'tra')
        % construct the channels from a linear combination of all magnetometer coils
        lf = sens.tra * lf;
      end
 
    otherwise
      ft_error('unsupported volume conductor model for MEG');
  end % switch type for MEG

elseif iseeg
  switch ft_headmodeltype(headmodel)

    case 'multisphere'
      % Based on the approximation of the potential due to a single dipole in
      % a multishell sphere by three dipoles in a homogeneous sphere, code
      % contributed by Punita Christopher. Note that this one should not get
      % confused with the MEG localspheres model.

      Nelec    = size(sens.elecpos, 1);
      Nspheres = length(headmodel.r);

      % the center of the spherical volume conduction model does not have
      % to be in the origin, therefore shift the spheres, the electrodes
      % and the dipole
      if isfield(headmodel, 'o')
        center = headmodel.o;
      else
        center = [0 0 0];
      end

      % sort the spheres from the smallest to the largest
      % furthermore, the radius should be one (?)
      [radii, indx] = sort(headmodel.r/max(headmodel.r));
      sigma = headmodel.cond(indx);
      r = (sens.elecpos-repmat(center, Nelec, 1))./max(headmodel.r);
      dippos = dippos./max(headmodel.r);

      if Ndipoles>1
        % loop over multiple dipoles
        lf = zeros(Nelec, 3*Ndipoles);
        for i=1:Ndipoles
          rq = dippos(i, :) - center;
          % compute the potential for each dipole ortientation
          % it would be much more efficient to change the punita function
          q1 = [1 0 0]; lf(:, (3*i-2)) = multisphere(Nspheres, radii, sigma, r, rq, q1);
          q1 = [0 1 0]; lf(:, (3*i-1)) = multisphere(Nspheres, radii, sigma, r, rq, q1);
          q1 = [0 0 1]; lf(:, (3*i )) = multisphere(Nspheres, radii, sigma, r, rq, q1);
        end
      else
        % only single dipole
        lf = zeros(Nelec, 3);
        rq = dippos - center;
        % compute the potential for each dipole ortientation
        % it would be much more efficient to change the punita function
        q1 = [1 0 0] ; lf(:, 1) = multisphere(Nspheres, radii, sigma, r, rq, q1);
        q1 = [0 1 0] ; lf(:, 2) = multisphere(Nspheres, radii, sigma, r, rq, q1);
        q1 = [0 0 1] ; lf(:, 3) = multisphere(Nspheres, radii, sigma, r, rq, q1);
      end

    case {'singlesphere', 'concentricspheres'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % EEG spherical volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % FIXME, this is not consistent between spherical and BEM
      % sort the spheres from the smallest to the largest
      [headmodel.r, indx] = sort(headmodel.r);
      headmodel.cond = headmodel.cond(indx);

      Nspheres = length(headmodel.cond);
      if length(headmodel.r)~=Nspheres
        ft_error('the number of spheres in the volume conductor model is ambiguous');
      end

      if isfield(headmodel, 'o')
        % shift the origin of the spheres, electrodes and dipole
        sens.elecpos = sens.elecpos - repmat(headmodel.o, size(sens.elecpos, 1), 1);
        dippos = dippos - repmat(headmodel.o, Ndipoles, 1);
      end

      switch Nspheres
        case 1
          funnam = 'eeg_leadfield1';
        case 2
          headmodel.r = [headmodel.r(1) headmodel.r(2) headmodel.r(2) headmodel.r(2)];
          headmodel.cond = [headmodel.cond(1) headmodel.cond(2) headmodel.cond(2) headmodel.cond(2)];
          funnam = 'eeg_leadfield4';
        case 3
          headmodel.r = [headmodel.r(1) headmodel.r(2) headmodel.r(3) headmodel.r(3)];
          headmodel.cond = [headmodel.cond(1) headmodel.cond(2) headmodel.cond(3) headmodel.cond(3)];
          funnam = 'eeg_leadfield4';
        case 4
          headmodel.r = [headmodel.r(1) headmodel.r(2) headmodel.r(3) headmodel.r(4)];
          headmodel.cond = [headmodel.cond(1) headmodel.cond(2) headmodel.cond(3) headmodel.cond(4)];
          funnam = 'eeg_leadfield4';
        otherwise
          ft_error('more than 4 concentric spheres are not supported')
      end

      lf = zeros(size(sens.elecpos, 1), 3*Ndipoles);
      for i=1:Ndipoles
        lf(:, (3*i-2):(3*i)) = feval(funnam, dippos(i, :), sens.elecpos, headmodel);
      end

    case {'bem', 'dipoli', 'asa', 'bemcp'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % EEG boundary element method volume conductor model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      lf = eeg_leadfieldb(dippos, sens.elecpos, headmodel);

    case 'openmeeg'
      ft_hastoolbox('openmeeg', 1);

      dsm         = ft_getopt(varargin, 'dsm');
      nonadaptive = ft_getopt(varargin, 'nonadaptive');

      [h2sens,ds2sens] = ft_sensinterp_openmeeg(dippos, headmodel, sens);
      if isempty(dsm)
        dsm            = ft_sysmat_openmeeg(dippos, headmodel, sens, nonadaptive);
      end
      lf               = ds2sens + h2sens*headmodel.mat*dsm;

    case {'infinite_currentdipole' 'infinite'}
      lf = eeg_infinite_dipole(dippos, sens.elecpos, headmodel);

    case 'halfspace'
      lf = eeg_halfspace_dipole(dippos, sens.elecpos, headmodel);

    case 'infinite_monopole'
      lf = eeg_infinite_monopole(dippos, sens.elecpos, headmodel);

    case 'halfspace_monopole'
      lf = eeg_halfspace_monopole(dippos, sens.elecpos, headmodel);

    case 'slab_monopole'
      lf = eeg_slab_monopole(dippos, sens.elecpos, headmodel);

    case 'simbio'
      ft_hastoolbox('simbio', 1);
      % note that the electrode information is contained in the headmodel (thanks to ft_prepare_vol_sens)
      lf = leadfield_simbio(dippos, headmodel);


    case 'duneuro'
      ft_hastoolbox('duneuro', 1);

      % note that the electrode information is contained in the headmodel
      lf = leadfield_duneuro(dippos, headmodel, 'eeg');

    case 'metufem'
      p3 = zeros(Ndipoles * 3, 6);
      for i = 1:Ndipoles
        p3((3*i - 2) : (3 * i), 1:3) = [dippos(i, :); dippos(i, :); dippos(i, :)];
        p3((3*i - 2) : (3 * i), 4:6) = [1 0 0; 0 1 0; 0 0 1];
      end
      lf = metufem('pot', p3', 'interp');

    case 'metubem'
      session = headmodel.session;
      p3 = zeros(Ndipoles * 3, 6);
      for i = 1:Ndipoles
        p3((3*i - 2) : (3 * i), 1:3) = [dippos(i, :); dippos(i, :); dippos(i, :)];
        p3((3*i - 2) : (3 * i), 4:6) = [1 0 0; 0 1 0; 0 0 1];
      end
      [lf, session] = bem_solve_lfm_eeg(session, p3);

    case 'fns'
      % note that the electrode information is contained in the headmodel
      % tolerance = 1e-8;
      lf = leadfield_fns(dippos, headmodel);

    case 'interpolate'
      % note that the electrode information is contained in the headmodel
      lf = leadfield_interpolate(dippos, headmodel);
      % the leadfield is already correctly referenced, i.e. it represents the
      % channel values rather than the electrode values. Prevent that the
      % referencing is done once more.
      sens.tra = speye(length(headmodel.filename));

    otherwise
      ft_error('unsupported volume conductor model for EEG');

  end % switch type for EEG

  % the forward model potential is computed on the electrodes relative to
  % an unknown reference, not on the channels. Therefore the data has to be
  % explicitly referenced here.
  if isfield(sens, 'tra')
    % apply the correct montage to the leadfield
    lf = sens.tra*lf;
  else
    % compute average reference for EEG leadfield
    for i=1:size(lf,2)
      lf(:,i) = lf(:,i) - mean(lf(:,i));
    end
  end

end % iseeg or ismeg

% optionally apply leadfield rank reduction
switch reducerank
  case 'yes'
    reducerank = 2;
  case 'no'
    reducerank = 3;
  otherwise
    % assume that it is specified as a number, keep it like this
end

if reducerank<size(lf,2)
  % decompose the leadfield
  for ii=1:Ndipoles
    tmplfd=lf(:, (3*ii-2):(3*ii));
    [u, s, v] = svd(tmplfd);
    r = diag(s);
    s(:) = 0;
    for j=1:reducerank
      s(j, j) = r(j);
    end

    if istrue(backproject)
      % recompose the leadfield with reduced rank
      lf(:, (3*ii-2):(3*ii)) = u * s * v';
    else
      % if not backprojected, the new leadfield has a different dimension
      if ii==1
        newlf    = zeros(size(lf,1), Ndipoles*reducerank);
        origrank = size(lf,2)./Ndipoles;
      end
      newlf(:, reducerank*(ii-1) + (1:reducerank)) = lf(:, origrank*(ii-1) + (1:origrank))*v(:,1:reducerank);
    end
  end

  if ~istrue(backproject)
    lf = newlf;
  end
  clear newlf;
end

% optionally apply leadfield normalization
switch normalize
  case 'yes'
    for ii=1:Ndipoles
      tmplf = lf(:, (3*ii-2):(3*ii));
      if normalizeparam==0.5
        % normalize the leadfield by the Frobenius norm of the matrix
        % this is the same as below in case normalizeparam is 0.5
        nrm = norm(tmplf, 'fro');
      else
        % normalize the leadfield by sum of squares of the elements of the leadfield matrix to the power "normalizeparam"
        % this is the same as the Frobenius norm if normalizeparam is 0.5
        nrm = sum(tmplf(:).^2)^normalizeparam;
      end
      if nrm>0
        tmplf = tmplf ./ nrm;
      end
      lf(:, (3*ii-2):(3*ii)) = tmplf;
    end
  case 'column'
    % normalize each column of the leadfield by its norm
    for ii=1:Ndipoles
      tmplf = lf(:, (3*ii-2):(3*ii));
      for j=1:size(tmplf, 2)
        nrm = sum(tmplf(:, j).^2)^normalizeparam;
        tmplf(:, j) = tmplf(:, j)./nrm;
      end
      lf(:, (3*ii-2):(3*ii)) = tmplf;
    end
end

% optionally apply a weight to the leadfield for each dipole location
if ~isempty(weight)
  for i=1:Ndipoles
    lf(:, 3*(i-1)+1) = lf(:, 3*(i-1)+1) * weight(i); % the leadfield for the x-direction
    lf(:, 3*(i-1)+2) = lf(:, 3*(i-1)+2) * weight(i); % the leadfield for the y-direction
    lf(:, 3*(i-1)+3) = lf(:, 3*(i-1)+3) * weight(i); % the leadfield for the z-direction
  end
end

if ~isempty(chanunit) || ~isempty(dipoleunit)
  assert(strcmp(headmodel.unit,  'm'), 'unit conversion only possible for SI input units');
  assert(strcmp(sens.unit, 'm'), 'unit conversion only possible for SI input units');
end

if ~isempty(chanunit)
  assert(all(strcmp(sens.chanunit, 'V') | strcmp(sens.chanunit, 'V/m') | strcmp(sens.chanunit, 'T') | strcmp(sens.chanunit, 'T/m')), 'unit conversion only possible for SI input units');
  % compute conversion factor and multiply each row of the matrix
  scale = cellfun(@ft_scalingfactor, sens.chanunit(:), chanunit(:));
  lf = bsxfun(@times, lf, scale(:));
  % prior to this conversion, the units might be  (T/m)/(A*m) for planar gradients or   (V/m)/(A*m) for bipolar EEG
  % after this conversion, the units will be     (T/cm)/(A*m)                      or (uV/mm)/(A*m)
end

if ~isempty(dipoleunit)
  scale = ft_scalingfactor('A*m', dipoleunit); % compue the scaling factor from A*m to the desired dipoleunit
  lf    = lf/scale;                         % the leadfield is expressed in chanunit per dipoleunit, i.e. chanunit/dipoleunit
end
