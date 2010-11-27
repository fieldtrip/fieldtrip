function [vol, sens] = ft_prepare_vol_sens(vol, sens, varargin)

% FT_PREPARE_VOL_SENS does some bookkeeping to ensure that the volume
% conductor model and the sensor array are ready for subsequent forward
% leadfield computations. It takes care of some pre-computations that can
% be done efficiently prior to the leadfield calculations.
%
% Use as
%   [vol, sens] = ft_prepare_vol_sens(vol, sens, ...)
% with input arguments
%   sens   structure with gradiometer or electrode definition
%   vol    structure with volume conductor definition
%
% The vol structure represents a volume conductor model, its contents
% depend on the type of model. The sens structure represents a sensor
% array, i.e. EEG electrodes or MEG gradiometers.
%
% Additional options should be specified in key-value pairs and can be
%   'channel'    cell-array with strings (default = 'all')
%   'order'      number, for single shell "Nolte" model (default = 10)
%
% The detailled behaviour of this function depends on whether the input
% consists of EEG or MEG and furthermoree depends on the type of volume
% conductor model:
% - in case of EEG single and concentric sphere models, the electrodes are
%   projected onto the skin surface.
% - in case of EEG boundary element models, the electrodes are projected on
%   the surface and a blilinear interpoaltion matrix from vertices to
%   electrodes is computed.
% - in case of MEG and a multispheres model, a local sphere is determined
%   for each coil in the gradiometer definition.
%  - in case of MEG with a singleshell Nolte model, the volume conduction
%    model is initialized
% In any case channel selection and reordering will be done. The channel
% order returned by this function corresponds to the order in the 'channel'
% option, or if not specified, to the order in the input sensor array.
%
% See also FT_READ_VOL, FT_READ_SENS, FT_TRANSFORM_VOL, FT_TRANSFORM_SENS, FT_COMPUTE_LEADFIELD

% Copyright (C) 2004-2009, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% get the options
% fileformat = keyval('fileformat',  varargin);
channel = keyval('channel',  varargin);  % cell-array with channel labels
order   = keyval('order',    varargin);  % order of expansion for Nolte method; 10 should be enough for real applications; in simulations it makes sense to go higher

% set the defaults
if isempty(channel),  channel = sens.label;   end
if isempty(order),    order = 10;             end

% determine whether the input contains EEG or MEG sensors
iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');

% determine the skin compartment
if ~isfield(vol, 'skin')
  if isfield(vol, 'bnd')
    vol.skin   = find_outermost_boundary(vol.bnd);
  elseif isfield(vol, 'r') && length(vol.r)<=4
    [dum, vol.skin] = max(vol.r);
  end
end

% determine the brain compartment
if ~isfield(vol, 'brain')
  if isfield(vol, 'bnd')
    vol.brain  = find_innermost_boundary(vol.bnd);
  elseif isfield(vol, 'r') && length(vol.r)<=4
    [dum, vol.brain] = min(vol.r);
  end
end

% otherwise the voltype assignment to an empty struct below won't work
if isempty(vol)
  vol = [];
end

% this makes them easier to recognise
sens.type = ft_senstype(sens);
vol.type  = ft_voltype(vol);


if isfield(vol, 'unit') && isfield(sens, 'unit') && ~strcmp(vol.unit, sens.unit)
  error('inconsistency in the units of the volume conductor and the sensor array');
end

if ismeg && iseeg
  % this is something that could be implemented relatively easily
  error('simultaneous EEG and MEG not yet supported');

elseif ~ismeg && ~iseeg
  error('the input does not look like EEG, nor like MEG');

elseif ismeg
  % always ensure that there is a linear transfer matrix for combining the coils into gradiometers
  if ~isfield(sens, 'tra');
    sens.tra = sparse(eye(length(sens.label)));
  end

  % select the desired channels from the gradiometer array
  % order them according to the users specification
  [selchan, selsens] = match_str(channel, sens.label);

  % first only modify the linear combination of coils into channels
  sens.label = sens.label(selsens);
  sens.tra   = sens.tra(selsens,:);
  % subsequently remove the coils that do not contribute to any sensor output
  selcoil  = find(sum(sens.tra,1)~=0);
  sens.pnt = sens.pnt(selcoil,:);
  sens.ori = sens.ori(selcoil,:);
  sens.tra = sens.tra(:,selcoil);

  switch ft_voltype(vol)
    case 'infinite'
      % nothing to do

    case 'singlesphere'
      % nothing to do

    case 'concentric'
      % nothing to do

    case 'neuromag'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if the forward model is computed using the external Neuromag toolbox,
      % we have to add a selection of the channels so that the channels
      % in the forward model correspond with those in the data.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [selchan, selsens] = match_str(channel, sens.label);
      vol.chansel = selsens;

    case 'multisphere'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % If the volume conduction model consists of multiple spheres then we
      % have to match the channels in the gradiometer array and the volume
      % conduction model.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % the initial multisphere volume conductor has a single sphere per
      % channel, whereas it should have a sphere for each coil
      if size(vol.r,1)==size(sens.pnt,1) && ~isfield(vol, 'label')
        % it appears that each coil already has a sphere, which suggests
        % that the volume conductor already has been prepared to match the
        % sensor array
        return
      elseif size(vol.r,1)==size(sens.pnt,1) && isfield(vol, 'label')
        if ~isequal(vol.label(:), sens.label(:))
          % if only the order is different, it would be possible to reorder them
          error('the coils in the volume conduction model do not correspond to the sensor array');
        else
          % the coil-specific spheres in the volume conductor should not have a label
          % because the label is already specified for the coils in the
          % sensor array
          vol = rmfield(vol, 'label');
        end
        return
      end

      % select the desired channels from the multisphere volume conductor
      % order them according to the users specification
      [selchan, selvol] = match_str(channel, vol.label);
      vol.label = vol.label(selvol);
      vol.r     = vol.r(selvol);
      vol.o     = vol.o(selvol,:);

      % the CTF way of representing the headmodel is one-sphere-per-channel
      % whereas the FieldTrip way of doing the forward computation is one-sphere-per-coil
      Nchans   = size(sens.tra,1);
      Ncoils   = size(sens.tra,2);
      Nspheres = size(vol.label);

      if isfield(vol, 'orig')
        % these are present in a CTF *.hdm file
        singlesphere.o(1,1) = vol.orig.MEG_Sphere.ORIGIN_X;
        singlesphere.o(1,2) = vol.orig.MEG_Sphere.ORIGIN_Y;
        singlesphere.o(1,3) = vol.orig.MEG_Sphere.ORIGIN_Z;
        singlesphere.r      = vol.orig.MEG_Sphere.RADIUS;
        % ensure consistent units
        singlesphere = ft_convert_units(singlesphere, vol.unit);
      else
        singlesphere = [];
      end

      if ~isempty(singlesphere)
        % determine the channels that do not have a corresponding sphere
        % and use the globally fitted single sphere for those
        missing = setdiff(sens.label, vol.label);
        if ~isempty(missing)
          warning('using the global fitted single sphere for %d channels that do not have a local sphere', length(missing));
        end
        for i=1:length(missing)
          vol.label(end+1) = missing(i);
          vol.r(end+1,:)   = singlesphere.r;
          vol.o(end+1,:)   = singlesphere.o;
        end
      end

      multisphere = [];
      % for each coil in the MEG helmet, determine the corresponding local sphere
      for i=1:Ncoils
        coilindex = find(sens.tra(:,i)~=0); % to which channel does the coil belong
        if length(coilindex)>1
          % this indicates that there are multiple channels to which this coil contributes,
          % which happens if the sensor array represents a synthetic higher-order gradient.
          [dum, coilindex] = max(abs(sens.tra(:,i)));
        end

        coillabel = sens.label{coilindex};  % what is the label of the channel
        chanindex = strmatch(coillabel, vol.label, 'exact');
        multisphere.r(i,:) = vol.r(chanindex);
        multisphere.o(i,:) = vol.o(chanindex,:);
      end
      vol = multisphere;

    case 'nolte'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if the forward model is computed using the code from Guido Nolte, we
      % have to initialize the volume model using the gradiometer coil
      % locations
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % compute the surface normals for each vertex point
      if ~isfield(vol.bnd, 'nrm')
        fprintf('computing surface normals\n');
        vol.bnd.nrm = normals(vol.bnd.pnt, vol.bnd.tri);
      end

      % estimate center and radius
      [center,radius] = fitsphere(vol.bnd.pnt);

      % initialize the forward calculation (only if gradiometer coils are available)
      if size(sens.pnt,1)>0
        vol.forwpar = meg_ini([vol.bnd.pnt vol.bnd.nrm], center', order, [sens.pnt sens.ori]);
      end

    case 'openmeeg'
        % nothing ?
        
    otherwise
      error('unsupported volume conductor model for MEG');
  end

elseif iseeg
  % select the desired channels from the electrode array
  % order them according to the users specification
  [selchan, selsens] = match_str(channel, sens.label);
  sens.label = sens.label(selsens);
  sens.pnt   = sens.pnt(selsens,:);

  % create a 2D projection and triangulation
  sens.prj   = elproj(sens.pnt);
  sens.tri   = delaunay(sens.prj(:,1), sens.prj(:,2));

  switch ft_voltype(vol)
    case 'infinite'
      % nothing to do

    case {'singlesphere', 'concentric'}
      % ensure that the electrodes ly on the skin surface
      radius = max(vol.r);
      pnt    = sens.pnt;
      if isfield(vol, 'o')
        % shift the the centre of the sphere to the origin
        pnt(:,1) = pnt(:,1) - vol.o(1);
        pnt(:,2) = pnt(:,2) - vol.o(2);
        pnt(:,3) = pnt(:,3) - vol.o(3);
      end
      distance = sqrt(sum(pnt.^2,2)); % to the center of the sphere
      if any((abs(distance-radius)/radius)>0.005)
        warning('electrodes do not lie on skin surface -> using radial projection')
      end
      pnt = pnt * radius ./ [distance distance distance];
      if isfield(vol, 'o')
        % shift the center back to the original location
        pnt(:,1) = pnt(:,1) + vol.o(1);
        pnt(:,2) = pnt(:,2) + vol.o(2);
        pnt(:,3) = pnt(:,3) + vol.o(3);
      end
      sens.pnt = pnt;

    case {'bem', 'dipoli', 'asa', 'avo', 'bemcp', 'openmeeg'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % do postprocessing of volume and electrodes in case of BEM model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % project the electrodes on the skin and determine the bilinear interpolation matrix
      if ~isfield(vol, 'tra')
        % determine boundary corresponding with skin and brain
        if ~isfield(vol, 'skin')
          vol.skin = find_outermost_boundary(vol.bnd);
          fprintf('determining skin compartment (%d)\n', vol.skin);
        end
        if ~isfield(vol, 'source')
          vol.source = find_innermost_boundary(vol.bnd);
          fprintf('determining source compartment (%d)\n', vol.source);
        end
        if size(vol.mat,1)~=size(vol.mat,2) && size(vol.mat,1)==length(sens.pnt)
          fprintf('electrode transfer and system matrix were already combined\n');
        else
          fprintf('projecting electrodes on skin surface\n');
          % compute linear interpolation from triangle vertices towards electrodes
          [el, prj] = project_elec(sens.pnt, vol.bnd(vol.skin).pnt, vol.bnd(vol.skin).tri);
          tra       = transfer_elec(vol.bnd(vol.skin).pnt, vol.bnd(vol.skin).tri, el);
          
          % replace the original electrode positions by the projected positions
          sens.pnt = prj;

          if size(vol.mat,1)==size(vol.bnd(vol.skin).pnt,1)
            % construct the transfer from only the skin vertices towards electrodes
            interp = tra;
          else
            % construct the transfer from all vertices (also brain/skull) towards electrodes
            interp = [];
            for i=1:length(vol.bnd)
              if i==vol.skin
                interp = [interp, tra];
              else
                interp = [interp, zeros(size(el,1), size(vol.bnd(i).pnt,1))];
              end
            end
          end

          % incorporate the linear interpolation matrix and the system matrix into one matrix
          % this speeds up the subsequent repeated leadfield computations
          fprintf('combining electrode transfer and system matrix\n');
          if strcmp(ft_voltype(vol), 'openmeeg')
            nb_points_external_surface = size(vol.bnd(vol.skin).pnt,1);
            vol.mat = vol.mat((end-nb_points_external_surface+1):end,:);            
            vol.mat = interp(:,1:nb_points_external_surface) * vol.mat;
          else
            % convert to sparse matrix to speed up the subsequent multiplication
            interp  = sparse(interp);
            vol.mat = interp * vol.mat;
            % ensure that the model potential will be average referenced
            avg = mean(vol.mat, 1);
            vol.mat = vol.mat - repmat(avg, size(vol.mat,1), 1);
          end
        end
      end
      
    otherwise
      error('unsupported volume conductor model for EEG');
  end

  % FIXME this needs carefull thought to ensure that the average referencing which is now done here and there, and that the linear interpolation in case of BEM are all dealt with consistently
  % % always ensure that there is a linear transfer matrix for
  % % rereferencing the EEG potential
  % if ~isfield(sens, 'tra');
  %   sens.tra = sparse(eye(length(sens.label)));
  % end

end % if iseeg or ismeg

