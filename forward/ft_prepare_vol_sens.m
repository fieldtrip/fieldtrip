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
% The detailed behaviour of this function depends on whether the input
% consists of EEG or MEG and furthermoree depends on the type of volume
% conductor model:
% - in case of EEG single and concentric sphere models, the electrodes are
%   projected onto the skin surface.
% - in case of EEG boundary element models, the electrodes are projected on
%   the surface and a blilinear interpoaltion matrix from vertices to
%   electrodes is computed.
% - in case of MEG and a localspheres model, a local sphere is determined
%   for each coil in the gradiometer definition.
%  - in case of MEG with a singleshell Nolte model, the volume conduction
%    model is initialized
% In any case channel selection and reordering will be done. The channel
% order returned by this function corresponds to the order in the 'channel'
% option, or if not specified, to the order in the input sensor array.
%
% See also FT_COMPUTE_LEADFIELD, FT_READ_VOL, FT_READ_SENS, FT_TRANSFORM_VOL,
% FT_TRANSFORM_SENS

% Copyright (C) 2004-2013, Robert Oostenveld
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

% get the optional input arguments
% fileformat = ft_getopt(varargin, 'fileformat');
channel = ft_getopt(varargin, 'channel', sens.label);   % cell-array with channel labels, default is all
order   = ft_getopt(varargin, 'order', 10);             % order of expansion for Nolte method; 10 should be enough for real applications; in simulations it makes sense to go higher

% ensure that the sensor description is up-to-date (Aug 2011)
sens = ft_datatype_sens(sens);

% this is to support volumes saved in mat-files, particularly interpolated
if ischar(vol)
  vpath = fileparts(vol);   % remember the path to the file
  vol   = ft_read_vol(vol); % replace the filename with the content of the file
end

% ensure that the volume conduction description is up-to-date (Jul 2012)
vol = ft_datatype_headmodel(vol);

% determine whether the input contains EEG or MEG sensors
iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');

% determine the skin compartment
if ~isfield(vol, 'skin_surface')
  if isfield(vol, 'bnd')
    vol.skin_surface   = find_outermost_boundary(vol.bnd);
  elseif isfield(vol, 'r') && length(vol.r)<=4
    [dum, vol.skin_surface] = max(vol.r);
  end
end

% determine the inner_skull_surface compartment
if ~isfield(vol, 'inner_skull_surface')
  if isfield(vol, 'bnd')
    vol.inner_skull_surface  = find_innermost_boundary(vol.bnd);
  elseif isfield(vol, 'r') && length(vol.r)<=4
    [dum, vol.inner_skull_surface] = min(vol.r);
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
  
  % keep a copy of the original sensor array, this is needed for the MEG localspheres model
  sens_orig = sens;
  
  % always ensure that there is a linear transfer matrix for combining the coils into gradiometers
  if ~isfield(sens, 'tra');
    Nchans = length(sens.label);
    Ncoils = size(sens.coilpos,1);
    if Nchans~=Ncoils
      error('inconsistent number of channels and coils');
    end
    sens.tra = eye(Nchans, Ncoils);
  end
  
  % select the desired channels from the gradiometer array
  % order them according to the users specification
  [selchan, selsens] = match_str(channel, sens.label);
  
  % first only modify the linear combination of coils into channels
  try, sens.chantype = sens.chantype(selsens,:); end
  try, sens.chanunit = sens.chanunit(selsens,:); end
  sens.chanpos  = sens.chanpos(selsens,:);
  sens.chanori  = sens.chanori(selsens,:);
  sens.label    = sens.label(selsens);
  sens.tra      = sens.tra(selsens,:);
  % subsequently remove the coils that do not contribute to any channel output
  selcoil      = any(sens.tra~=0,1);
  sens.coilpos = sens.coilpos(selcoil,:);
  sens.coilori = sens.coilori(selcoil,:);
  sens.tra     = sens.tra(:,selcoil);
  
  switch ft_voltype(vol)
    case {'infinite' 'infinite_monopole'}
      % nothing to do
      
    case 'singlesphere'
      % nothing to do
      
    case 'concentricspheres'
      % nothing to do
      
    case 'neuromag'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if the forward model is computed using the external Neuromag toolbox,
      % we have to add a selection of the channels so that the channels
      % in the forward model correspond with those in the data.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [selchan, selsens] = match_str(channel, sens.label);
      vol.chansel = selsens;
      
    case 'localspheres'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % If the volume conduction model consists of multiple spheres then we
      % have to match the channels in the gradiometer array and the volume
      % conduction model.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % use the original sensor array instead of the one with a subset of
      % channels, because we need the complete mapping of coils to channels
      sens = sens_orig;
      
      % remove the coils that do not contribute to any channel output
      % since these do not have a corresponding sphere
      selcoil      = find(sum(sens.tra,1)~=0);
      sens.coilpos = sens.coilpos(selcoil,:);
      sens.coilori = sens.coilori(selcoil,:);
      sens.tra     = sens.tra(:,selcoil);
      
      % the initial localspheres volume conductor has a local sphere per
      % channel, whereas it should have a local sphere for each coil
      if size(vol.r,1)==size(sens.coilpos,1) && ~isfield(vol, 'label')
        % it appears that each coil already has a sphere, which suggests
        % that the volume conductor already has been prepared to match the
        % sensor array
        return
      elseif size(vol.r,1)==size(sens.coilpos,1) && isfield(vol, 'label')
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
      
      localspheres = [];
      % for each coil in the MEG helmet, determine the corresponding channel and from that the corresponding local sphere
      for i=1:Ncoils
        coilindex = find(sens.tra(:,i)~=0); % to which channel does this coil belong
        if length(coilindex)>1
          % this indicates that there are multiple channels to which this coil contributes,
          % which happens if the sensor array represents a synthetic higher-order gradient.
          [dum, coilindex] = max(abs(sens.tra(:,i)));
        end
        
        coillabel = sens.label{coilindex};               % what is the label of this channel
        chanindex = find(strcmp(coillabel, vol.label));  % what is the index of this channel in the list of local spheres
        localspheres.r(i,:) = vol.r(chanindex);
        localspheres.o(i,:) = vol.o(chanindex,:);
      end
      vol = localspheres;
      
      % finally do the selection of channels and coils
      % order them according to the users specification
      [selchan, selsens] = match_str(channel, sens.label);
      
      % first only modify the linear combination of coils into channels
      try, sens.chantype = sens.chantype(selsens,:); end
      try, sens.chanunit = sens.chanunit(selsens,:); end
      sens.chanpos = sens.chanpos(selsens,:);
      sens.chanori = sens.chanori(selsens,:);
      sens.label   = sens.label(selsens);
      sens.tra     = sens.tra(selsens,:);
      % subsequently remove the coils that do not contribute to any sensor output
      selcoil      = find(sum(sens.tra,1)~=0);
      sens.coilpos = sens.coilpos(selcoil,:);
      sens.coilori = sens.coilori(selcoil,:);
      sens.tra     = sens.tra(:,selcoil);
      % make the same selection of coils in the localspheres model
      vol.r = vol.r(selcoil);
      vol.o = vol.o(selcoil,:);
      
    case 'singleshell'
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
      if size(sens.coilpos,1)>0 && ~isfield(vol, 'forwpar')
        vol.forwpar = meg_ini([vol.bnd.pnt vol.bnd.nrm], center', order, [sens.coilpos sens.coilori]);
      end
      
    case 'openmeeg'
      error('MEG not yet supported with openmeeg');
      
    case 'simbio'
      error('MEG not yet supported with simbio');
      
    otherwise
      error('unsupported volume conductor model for MEG');
  end
  
elseif iseeg
  
  % select the desired channels from the electrode array
  % order them according to the users specification
  [selchan, selsens] = match_str(channel, sens.label);
  Nchans = length(sens.label);
  
  if isfield(sens, 'tra')
    % first only modify the linear combination of electrodes into channels
    sens.chanpos = sens.chanpos(selsens,:);
    sens.label   = sens.label(selsens);
    sens.tra     = sens.tra(selsens,:);
    % subsequently remove the electrodes that do not contribute to any channel output
    selelec      = any(sens.tra~=0,1);
    sens.elecpos = sens.elecpos(selelec,:);
    sens.tra     = sens.tra(:,selelec);
  else
    % the electrodes and channels are identical
    sens.chanpos = sens.chanpos(selsens,:);
    sens.elecpos = sens.elecpos(selsens,:);
    sens.label   = sens.label(selsens);
  end
  
  switch ft_voltype(vol)
    case {'infinite' 'infinite_monopole'}
      % nothing to do
      
    case {'halfspace', 'halfspace_monopole'}
      % electrodes' all-to-all distances
      numelec = size(sens.elecpos,1);
      ref_el = sens.elecpos(1,:);
      md = dist( (sens.elecpos-repmat(ref_el,[numelec 1]))' );
      % take the min distance as reference
      md = min(md(1,2:end));
      pnt = sens.elecpos;
      % scan the electrodes and reposition the ones which are in the
      % wrong halfspace (projected on the plane)... if not too far away!
      for i=1:size(pnt,1)
        P = pnt(i,:);
        is_in_empty = acos(dot(vol.ori,(P-vol.pnt)./norm(P-vol.pnt))) < pi/2;
        if is_in_empty
          dPplane = abs(dot(vol.ori, vol.pnt-P, 2));
          if dPplane>md
            error('Some electrodes are too distant from the plane: consider repositioning them')
          else
            % project point on plane
            Ppr = pointproj(P,[vol.pnt vol.ori]);
            pnt(i,:) = Ppr;
          end
        end
      end
      sens.elecpos = pnt;
      
    case {'slab_monopole'}
      % electrodes' all-to-all distances
      numel  = size(sens.elecpos,1);
      ref_el = sens.elecpos(1,:);
      md  = dist( (sens.elecpos-repmat(ref_el,[numel 1]))' );
      % choose min distance between electrodes
      md  = min(md(1,2:end));
      pnt = sens.elecpos;
      % looks for contacts outside the strip which are not too far away
      % and projects them on the nearest plane
      for i=1:size(pnt,1)
        P = pnt(i,:);
        instrip1 = acos(dot(vol.ori1,(P-vol.pnt1)./norm(P-vol.pnt1))) > pi/2;
        instrip2 = acos(dot(vol.ori2,(P-vol.pnt2)./norm(P-vol.pnt2))) > pi/2;
        is_in_empty = ~(instrip1&instrip2);
        if is_in_empty
          dPplane1 = abs(dot(vol.ori1, vol.pnt1-P, 2));
          dPplane2 = abs(dot(vol.ori2, vol.pnt2-P, 2));
          if dPplane1>md && dPplane2>md
            error('Some electrodes are too distant from the planes: consider repositioning them')
          elseif dPplane2>dPplane1
            % project point on nearest plane
            Ppr = pointproj(P,[vol.pnt1 vol.ori1]);
            pnt(i,:) = Ppr;
          else
            % project point on nearest plane
            Ppr = pointproj(P,[vol.pnt2 vol.ori2]);
            pnt(i,:) = Ppr;
          end
        end
      end
      sens.elecpos = pnt;
      
    case {'singlesphere', 'concentricspheres'}
      % ensure that the electrodes ly on the skin surface
      radius = max(vol.r);
      pnt    = sens.elecpos;
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
      sens.elecpos = pnt;
      
    case {'bem', 'dipoli', 'asa', 'bemcp', 'openmeeg'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % do postprocessing of volume and electrodes in case of BEM model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % project the electrodes on the skin and determine the bilinear interpolation matrix
      if ~isfield(vol, 'tra')
        % determine boundary corresponding with skin and inner_skull_surface
        if ~isfield(vol, 'skin_surface')
          vol.skin_surface = find_outermost_boundary(vol.bnd);
          fprintf('determining skin compartment (%d)\n', vol.skin_surface);
        end
        if ~isfield(vol, 'source')
          vol.source = find_innermost_boundary(vol.bnd);
          fprintf('determining source compartment (%d)\n', vol.source);
        end
        if size(vol.mat,1)~=size(vol.mat,2) && size(vol.mat,1)==length(sens.elecpos)
          fprintf('electrode transfer and system matrix were already combined\n');
        else
          fprintf('projecting electrodes on skin surface\n');
          % compute linear interpolation from triangle vertices towards electrodes
          [el, prj] = project_elec(sens.elecpos, vol.bnd(vol.skin_surface).pnt, vol.bnd(vol.skin_surface).tri);
          tra       = transfer_elec(vol.bnd(vol.skin_surface).pnt, vol.bnd(vol.skin_surface).tri, el);
          
          % replace the original electrode positions by the projected positions
          sens.elecpos = prj;
          
          if size(vol.mat,1)==size(vol.bnd(vol.skin_surface).pnt,1)
            % construct the transfer from only the skin vertices towards electrodes
            interp = tra;
          else
            % construct the transfer from all vertices (also inner_skull_surface/outer_skull_surface) towards electrodes
            interp = [];
            for i=1:length(vol.bnd)
              if i==vol.skin_surface
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
            % check that the external toolbox is present
            ft_hastoolbox('openmeeg', 1);
            
            nb_points_external_surface = size(vol.bnd(vol.skin_surface).pnt,1);
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
      
    case 'fns'
      if isfield(vol,'bnd')
        [el, prj] = project_elec(sens.elecpos, vol.bnd.pnt, vol.bnd.tri);
        sens.tra = transfer_elec(vol.bnd.pnt, vol.bnd.tri, el);
        % replace the original electrode positions by the projected positions
        sens.elecpos = prj;
      end
      
    case 'simbio'
      % check that the external toolbox is present
      ft_hastoolbox('simbio', 1);
      
      % extract the outer surface
      bnd = mesh2edge(vol);
      for j=1:length(sens.label)
        d = bsxfun(@minus, bnd.pnt, sens.elecpos(j,:));
        [d, i] = min(sum(d.^2, 2));
        % replace the position of each electrode by the closest vertex
        sens.elecpos(j,:) = bnd.pnt(i,:);
      end
      
      if isfield(sens, 'chanpos')
        % this is invalid after the projection to the surface
        sens = rmfield(sens, 'chanpos');
      end
      
      vol.transfer = sb_transfer(vol,sens);
      
    case 'interpolate'
      % this is to allow moving leadfield files
      if ~exist(vol.filename{1}, 'file')
         for i = 1:length(vol.filename)
             [p, f, x] = fileparts(vol.filename{i});
             vol.filename{i} = fullfile(vpath, [f x]);
         end
      end
       
      if ~isfield(sens, 'tra') && isequal(sens.chanpos, sens.elecpos)
        sens.tra = eye(size(sens.chanpos,1));
      end
      
      if ~isfield(vol.sens, 'tra') && isequal(vol.sens.chanpos, vol.sens.elecpos)
        vol.sens.tra = eye(size(vol.sens.chanpos,1));
      end
      
      % the channel positions can be nan, for example for a bipolar montage
      match = isequal(sens.label, vol.sens.label)    & ...
        isequalwithequalnans(sens.tra, vol.sens.tra) & ...
        isequal(sens.elecpos, vol.sens.elecpos)      & ...
        isequalwithequalnans(sens.chanpos, vol.sens.chanpos);
      
      if match
        % the input sensor array matches precisely with the forward model
        % no further interpolation is needed
      else
        % interpolate the channels in the forward model to the desired channels
        filename = tempname;
        vol  = ft_headmodel_interpolate(filename, sens, vol);
        % update the sensor array with the one from the volume conductor
        sens = vol.sens;
      end % if recomputing interpolation
      
      % for the leadfield computations the @nifti object is used to map the image data into memory
      ft_hastoolbox('spm8up', 1);
      for i=1:length(vol.sens.label)
        % map each of the leadfield files into memory
        vol.chan{i} = nifti(vol.filename{i});
      end
      
    otherwise
      error('unsupported volume conductor model for EEG');
  end
  
  % FIXME this needs careful thought to ensure that the average referencing which is now done here and there, and that the linear interpolation in case of BEM are all dealt with consistently
  % % always ensure that there is a linear transfer matrix for
  % % rereferencing the EEG potential
  % if ~isfield(sens, 'tra');
  %   sens.tra = eye(length(sens.label));
  % end
  
end % if iseeg or ismeg

if isfield(sens, 'tra')
  if issparse(sens.tra) && size(sens.tra, 1)==1
    % this multiplication would result in a sparse leadfield, which is not what we want
    % the effect can be demonstrated as sparse(1)*rand(1,10), see also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1169#c7
    sens.tra = full(sens.tra);
  elseif ~issparse(sens.tra) && size(sens.tra, 1)>1
    % the multiplication of the "sensor" leadfield (electrode or coil) with the tra matrix to get the "channel" leadfield
    % is faster for most cases if the pre-multiplying weighting matrix is made sparse
    sens.tra = sparse(sens.tra);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ppr = pointproj(P,plane)
% projects a point on a plane
% plane(1:3) is a point on the plane
% plane(4:6) is the ori of the plane
Ppr  = [];
ori  = plane(4:6);
line = [P ori];
% get indices of line and plane which are parallel
par = abs(dot(plane(4:6), line(:,4:6), 2))<1e-14;
% difference between origins of plane and line
dp = plane(1:3) - line(:, 1:3);
% Divide only for non parallel vectors (DL)
t = dot(ori(~par,:), dp(~par,:), 2)./dot(ori(~par,:), line(~par,4:6), 2);
% compute coord of intersection point
Ppr(~par, :) = line(~par,1:3) + repmat(t,1,3).*line(~par,4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as a replacement for the dist function in the Neural
% Networks toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d] = dist(x)
n = size(x,2);
d = zeros(n,n);
for i=1:n
  for j=(i+1):n
    d(i,j) = sqrt(sum((x(:,i)-x(:,j)).^2));
    d(j,i) = d(i,j);
  end
end
