function [vol, sens] = prepare_vol_sens(vol, sens, varargin)

% PREPARE_VOL_SENS does some bookkeeping to ensure that the volume
% conductor model and the sensor array are ready for subsequent forward
% leadfield computations. It takes care of some pre-computations that can
% be done efficiently prior to the leadfield calculations.
%
% Use as
%   [vol, sens] = prepare_vol_sens(vol, sens, ...)
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
% See also READ_VOL, READ_SENS, TRANSFORM_VOL, TRANSFORM_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2004-2009, Robert Oostenveld
%
% $Log: prepare_vol_sens.m,v $
% Revision 1.21  2009/09/27 19:12:38  crimic
% 	wrapper adapted for openmeeg forward solution
%
% Revision 1.20  2009/09/21 11:12:51  roboos
% added openmeeg as supported voltype, thanks to Cristiano
%
% Revision 1.19  2009/05/29 11:50:34  vlalit
% Fixed a bug with wrong kind of brackets
%
% Revision 1.18  2009/05/25 11:50:40  roboos
% consistent handling of multiple spheres in case of ctf localspheres.hdm and fieldtrip prepare_localspheres, also in case of synthetic gradients.
%
% Revision 1.17  2009/05/25 08:06:37  roboos
% don't assign the channel labels to each coil for multisphere
% handle the input situation of an already prepared vol and sens (for multisphere)
%
% Revision 1.16  2009/05/18 15:57:23  roboos
% added the label of each coil to the multisphere model output
%
% Revision 1.15  2009/04/01 12:36:52  roboos
% use Taubin's method for fitting the sphere instead of Guido's iterative sphfit function
%
% Revision 1.14  2009/03/26 16:44:03  roboos
% allow 3rd order gradients iduring the construction of the localspheres model, requires that the hdm file contains a global sphere
%
% Revision 1.13  2009/03/23 21:15:18  roboos
% fixed bug for empty vol (inf medium magnetic dipole)
%
% Revision 1.12  2009/03/11 11:29:26  roboos
% ensure that the channel order in the sens and in the vol is consistent with  the user-specified channel-keyval argument
%
% Revision 1.11  2009/02/02 13:06:40  roboos
% added bemcp
% changed handling of vertex->electrode interpolation, now also possible if bem system matrix only describes the skin surface
%
% Revision 1.10  2009/01/19 12:13:49  roboos
% added code at the end to determine the brain and the skin compartment
%
% Revision 1.9  2008/12/24 10:34:15  roboos
% added two fprintf statements
%
% Revision 1.8  2008/09/29 09:56:04  release
% some spelling fixes, thanks to Karl
%
% Revision 1.7  2008/07/22 10:17:15  roboos
% replaced identical with strcmp
%
% Revision 1.6  2008/07/21 20:28:44  roboos
% added check on units (mm/cm/m) of the sensor array and volume conductor, give error if inconsistent
%
% Revision 1.5  2008/04/30 13:47:59  roboos
% project electrodes on scalp surface if sphere
% always ensure that grad.tra exists (not yet for elec)
%
% Revision 1.4  2008/04/15 20:36:21  roboos
% added explicit handling of various BEM implementations, i.e. for all voltype variants
%
% Revision 1.3  2008/04/10 11:00:29  roboos
% fixed some small but obvious bugs
%
% Revision 1.2  2008/04/09 20:37:32  roboos
% copied code over from ft version, not yet tested
%
% Revision 1.1  2008/03/06 09:30:36  roboos
% Created skeleton implementation according to how it should be for the forwinv toolbox, i.e. fieldtrip independent, so that it can be included in spm8.
% The functionality should be moved from the existing fieldtrip/private/prepare_vol_sens.m function into this new function.
%

% get the options
% fileformat = keyval('fileformat',  varargin);
channel = keyval('channel',  varargin);  % cell-array with channel labels
order   = keyval('order',    varargin);  % order of expansion for Nolte method; 10 should be enough for real applications; in simulations it makes sense to go higher

% set the defaults
if isempty(channel),  channel = sens.label;   end
if isempty(order),    order = 10;             end

% determine whether the input contains EEG or MEG sensors
iseeg = senstype(sens, 'eeg');
ismeg = senstype(sens, 'meg');

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

  switch voltype(vol)
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
      if size(vol.r,1)==size(sens.pnt,1)
        % it appears that each coil already has a sphere, which suggests
        % that the volume conductor already has been prepared to match the
        % sensor array
        if ~isfield(vol, 'label')
          % this is ok, since coils should not have labels
        else
          % hmm, this only works if there is a one-to-one match between
          % coils in the sensor array and coils in the volme conductor
          if ~isequal(vol.label(:), sens.label(:))
            % the problem here is that each channel can have multiple coils
            % in case of a gradiometer arrangement. The consequence is that
            % the coil label is not unique, because the bottom and top
            % coil will have the same label. That causes problems with the
            % channel selection.
            error('the coils in the volume conduction model do not correspond to the sensor array');
          else
            % the coil-specific spheres in the volume conductor should not have a label
            % because the label is already specified for the coils in the
            % sensor array
            vol = rmfield(vol, 'label');
          end
        end
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
        singlesphere = convert_units(singlesphere, vol.unit);
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

  switch voltype(vol)
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

    case {'bem', 'dipoli', 'asa', 'avo', 'bemcp'}
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
          el   = project_elec(sens.pnt, vol.bnd(vol.skin).pnt, vol.bnd(vol.skin).tri);
          tra  = transfer_elec(vol.bnd(vol.skin).pnt, vol.bnd(vol.skin).tri, el);
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
          % convert to sparse matrix to speed up the subsequent multiplication
          interp    = sparse(interp);
          % incorporate the linear interpolation matrix and the system matrix into one matrix
          % this speeds up the subsequent repeated leadfield computations
          fprintf('combining electrode transfer and system matrix\n');
          vol.mat = interp * vol.mat;
          % FIXME should I also add the electrode labels to the volume definition?
        end
        % ensure that the model potential will be average referenced
        avg = mean(vol.mat, 1);
        vol.mat = vol.mat - repmat(avg, size(vol.mat,1), 1);
      end

    case {'openmeeg'}
      % do nothing
      % in case of openmeeg do nothing because electrodes projection is
      % already performed in command line INRIA routines
      % FIXME: to be checked the average referencing of the openmeeg tool 
      % vol.mat = vol.mat;
      
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
sens.type = senstype(sens);
vol.type  = voltype(vol);
