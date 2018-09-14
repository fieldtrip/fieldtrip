function [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, varargin)

% FT_PREPARE_VOL_SENS does some bookkeeping to ensure that the volume
% conductor model and the sensor array are ready for subsequent forward
% leadfield computations. It takes care of some pre-computations that can
% be done efficiently prior to the leadfield calculations.
%
% Use as
%   [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, ...)
% with input arguments
%   headmodel  structure with volume conductor definition
%   sens       structure with gradiometer or electrode definition
%
% The headmodel structure represents a volume conductor model of the head,
% its contents depend on the type of model. The sens structure represents a
% sensor array, i.e. EEG electrodes or MEG gradiometers.
%
% Additional options should be specified in key-value pairs and can be
%   'channel'    cell-array with strings (default = 'all')
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

% Copyright (C) 2004-2015, Robert Oostenveld
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

if iscell(headmodel) && iscell(sens)
  % this represents combined EEG, ECoG and/or MEG
  for i=1:numel(headmodel)
    [headmodel{i}, sens{i}] = ft_prepare_vol_sens(headmodel{i}, sens{i}, varargin{:});
  end
  return
end

% get the optional input arguments
% fileformat = ft_getopt(varargin, 'fileformat');
channel = ft_getopt(varargin, 'channel', sens.label);   % cell-array with channel labels, default is all

% ensure that the sensor description is up-to-date (Aug 2011)
sens = ft_datatype_sens(sens);

% this is to support volumes saved in mat-files, particularly interpolated
if ischar(headmodel)
  vpath     = fileparts(headmodel);   % remember the path to the file
  headmodel = ft_read_vol(headmodel); % replace the filename with the content of the file
end

% ensure that the volume conduction description is up-to-date (Jul 2012)
headmodel = ft_datatype_headmodel(headmodel);

% determine whether the input contains EEG or MEG sensors
iseeg = ft_senstype(sens, 'eeg');
ismeg = ft_senstype(sens, 'meg');

% determine the skin compartment
if ~isfield(headmodel, 'skin_surface')
  if isfield(headmodel, 'bnd')
    headmodel.skin_surface   = find_outermost_boundary(headmodel.bnd);
  elseif isfield(headmodel, 'r') && length(headmodel.r)<=4
    [dum, headmodel.skin_surface] = max(headmodel.r);
  end
end

% determine the inner_skull_surface compartment
if ~isfield(headmodel, 'inner_skull_surface')
  if isfield(headmodel, 'bnd')
    headmodel.inner_skull_surface  = find_innermost_boundary(headmodel.bnd);
  elseif isfield(headmodel, 'r') && length(headmodel.r)<=4
    [dum, headmodel.inner_skull_surface] = min(headmodel.r);
  end
end

% otherwise the voltype assignment to an empty struct below won't work
if isempty(headmodel)
  headmodel = [];
end

% this makes them easier to recognise
sens.type = ft_senstype(sens);
headmodel.type  = ft_voltype(headmodel);

if isfield(headmodel, 'unit') && isfield(sens, 'unit') && ~strcmp(headmodel.unit, sens.unit)
  ft_error('inconsistency in the units of the volume conductor and the sensor array');
end

if ismeg && iseeg
  % this is something that could be implemented relatively easily
  ft_error('simultaneous EEG and MEG not yet supported');
  
elseif ~ismeg && ~iseeg
  ft_error('the input does not look like EEG, nor like MEG');
  
elseif ismeg
  
  % always ensure that there is a linear transfer matrix for combining the coils into gradiometers
  if ~isfield(sens, 'tra');
    Nchans = length(sens.label);
    Ncoils = size(sens.coilpos,1);
    if Nchans~=Ncoils
      ft_error('inconsistent number of channels and coils');
    end
    sens.tra = eye(Nchans, Ncoils);
  end
  
  if ~ft_voltype(headmodel, 'localspheres')
    % select the desired channels from the gradiometer array
    [selchan, selsens] = match_str(channel, sens.label);
    % only keep the desired channels, order them according to the users specification
    try, sens.chantype = sens.chantype(selsens,:); end
    try, sens.chanunit = sens.chanunit(selsens,:); end
    try, sens.chanpos  = sens.chanpos (selsens,:); end
    try, sens.chanori  = sens.chanori (selsens,:); end
    sens.label    = sens.label(selsens);
    sens.tra      = sens.tra(selsens,:);
  else
    % for the localspheres model it is done further down
  end
  
  % remove the coils that do not contribute to any channel output
  selcoil      = any(sens.tra~=0,1);
  sens.coilpos = sens.coilpos(selcoil,:);
  sens.coilori = sens.coilori(selcoil,:);
  sens.tra     = sens.tra(:,selcoil);
  
  switch ft_voltype(headmodel)
    case {'infinite' 'infinite_monopole' 'infinite_currentdipole' 'infinite_magneticdipole'}
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
      headmodel.chansel = selsens;
      
    case 'localspheres'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % If the volume conduction model consists of multiple spheres then we
      % have to match the channels in the gradiometer array and the volume
      % conduction model.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % the initial localspheres volume conductor has a local sphere per
      % channel, whereas it should have a local sphere for each coil
      if size(headmodel.r,1)==size(sens.coilpos,1) && ~isfield(headmodel, 'label')
        % it appears that each coil already has a sphere, which suggests
        % that the volume conductor already has been prepared to match the
        % sensor array
        return
      elseif size(headmodel.r,1)==size(sens.coilpos,1) && isfield(headmodel, 'label')
        if ~isequal(headmodel.label(:), sens.label(:))
          % if only the order is different, it would be possible to reorder them
          ft_error('the coils in the volume conduction model do not correspond to the sensor array');
        else
          % the coil-specific spheres in the volume conductor should not have a label
          % because the label is already specified for the coils in the
          % sensor array
          headmodel = rmfield(headmodel, 'label');
        end
        return
      end
      
      % the CTF way of representing the headmodel is one-sphere-per-channel
      % whereas the FieldTrip way of doing the forward computation is one-sphere-per-coil
      Nchans   = size(sens.tra,1);
      Ncoils   = size(sens.tra,2);
      Nspheres = size(headmodel.label);
      
      if isfield(headmodel, 'orig')
        % these are present in a CTF *.hdm file
        singlesphere.o(1,1) = headmodel.orig.MEG_Sphere.ORIGIN_X;
        singlesphere.o(1,2) = headmodel.orig.MEG_Sphere.ORIGIN_Y;
        singlesphere.o(1,3) = headmodel.orig.MEG_Sphere.ORIGIN_Z;
        singlesphere.r      = headmodel.orig.MEG_Sphere.RADIUS;
        % ensure consistent units
        singlesphere = ft_convert_units(singlesphere, headmodel.unit);
        % determine the channels that do not have a corresponding sphere
        % and use the globally fitted single sphere for those
        missing = setdiff(sens.label, headmodel.label);
        if ~isempty(missing)
          ft_warning('using the global fitted single sphere for %d channels that do not have a local sphere', length(missing));
        end
        for i=1:length(missing)
          headmodel.label(end+1) = missing(i);
          headmodel.r(end+1,:)   = singlesphere.r;
          headmodel.o(end+1,:)   = singlesphere.o;
        end
      end
      
      % make a new structure that only holds the local spheres, one per coil
      localspheres = [];
      localspheres.type = headmodel.type;
      localspheres.unit = headmodel.unit;
      
      % for each coil in the MEG helmet, determine the corresponding channel and from that the corresponding local sphere
      for i=1:Ncoils
        coilindex = find(sens.tra(:,i)~=0); % to which channel does this coil belong
        if length(coilindex)>1
          % this indicates that there are multiple channels to which this coil contributes,
          % which happens if the sensor array represents a synthetic higher-order gradient.
          [dum, coilindex] = max(abs(sens.tra(:,i)));
        end
        
        coillabel = sens.label{coilindex};               % what is the label of this channel
        chanindex = find(strcmp(coillabel, headmodel.label));  % what is the index of this channel in the list of local spheres
        localspheres.r(i,:) = headmodel.r(chanindex);
        localspheres.o(i,:) = headmodel.o(chanindex,:);
      end
      headmodel = localspheres;
      
      % finally do the selection of channels and coils
      % order them according to the users specification
      [selchan, selsens] = match_str(channel, sens.label);
      
      % first only modify the linear combination of coils into channels
      try, sens.chantype = sens.chantype(selsens,:); end
      try, sens.chanunit = sens.chanunit(selsens,:); end
      try, sens.chanpos  = sens.chanpos (selsens,:); end
      try, sens.chanori  = sens.chanori (selsens,:); end
      sens.label   = sens.label(selsens);
      sens.tra     = sens.tra(selsens,:);
      % subsequently remove the coils that do not contribute to any sensor output
      selcoil      = find(sum(sens.tra,1)~=0);
      sens.coilpos = sens.coilpos(selcoil,:);
      sens.coilori = sens.coilori(selcoil,:);
      sens.tra     = sens.tra(:,selcoil);
      % make the same selection of coils in the localspheres model
      headmodel.r = headmodel.r(selcoil);
      headmodel.o = headmodel.o(selcoil,:);
      
    case 'singleshell'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if the forward model is computed using the code from Guido Nolte, we
      % have to initialize the volume model using the gradiometer coil
      % locations
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % compute the surface normals for each vertex point
      if ~isfield(headmodel.bnd, 'nrm')
        fprintf('computing surface normals\n');
        headmodel.bnd.nrm = normals(headmodel.bnd.pos, headmodel.bnd.tri);
      end
      
      % estimate center and radius
      [center, radius] = fitsphere(headmodel.bnd.pos);

      % order of spherical spherical harmonics, for 'real' realistic volume conductors order=10 is o.k
      if isfield(headmodel, 'order')
        order = headmodel.order;
      else
        order = 10;
      end
      
      % initialize the forward calculation (only if  coils are available)
      if size(sens.coilpos,1)>0 && ~isfield(headmodel, 'forwpar')
        s = ft_scalingfactor(headmodel.unit, 'cm');
        headmodel.forwpar = meg_ini([s*headmodel.bnd.pos headmodel.bnd.nrm], s*center', order, [s*sens.coilpos sens.coilori]);
        headmodel.forwpar.scale = s;
      end
      
    case  'openmeeg' 
        % don't do anything, h2em or h2mm generated later in ft_prepare_leadfield

    case 'simbio'
      ft_error('MEG not yet supported with simbio');
      
    otherwise
      ft_error('unsupported volume conductor model for MEG');
  end
  
elseif iseeg
  
  % the electrodes are used, the channel positions are not relevant any more
  % channel positinos need to be recomputed after projecting the electrodes on the skin
  if isfield(sens, 'chanpos'); sens = rmfield(sens, 'chanpos'); end
  
  % select the desired channels from the electrode array
  % order them according to the users specification
  [selchan, selsens] = match_str(channel, sens.label);
  Nchans = length(sens.label);
  
  sens.label     = sens.label(selsens);
  try, sens.chantype  = sens.chantype(selsens); end
  try, sens.chanunit  = sens.chanunit(selsens); end
  
  if isfield(sens, 'tra')
    % first only modify the linear combination of electrodes into channels
    sens.tra     = sens.tra(selsens,:);
    % subsequently remove the electrodes that do not contribute to any channel output
    selelec      = any(sens.tra~=0,1);
    sens.elecpos = sens.elecpos(selelec,:);
    sens.tra     = sens.tra(:,selelec);
  else
    % the electrodes and channels are identical
    sens.elecpos = sens.elecpos(selsens,:);
  end
  
  switch ft_voltype(headmodel)
    case {'infinite' 'infinite_monopole' 'infinite_currentdipole'}
      % nothing to do
      
    case {'halfspace', 'halfspace_monopole'}
      % electrodes' all-to-all distances
      numelec = size(sens.elecpos,1);
      ref_el = sens.elecpos(1,:);
      md = dist( (sens.elecpos-repmat(ref_el,[numelec 1]))' );
      % take the min distance as reference
      md = min(md(1,2:end));
      pos = sens.elecpos;
      % scan the electrodes and reposition the ones which are in the
      % wrong halfspace (projected on the plane)... if not too far away!
      for i=1:size(pos,1)
        P = pos(i,:);
        is_in_empty = acos(dot(headmodel.ori,(P-headmodel.pos)./norm(P-headmodel.pos))) < pi/2;
        if is_in_empty
          dPplane = abs(dot(headmodel.ori, headmodel.pos-P, 2));
          if dPplane>md
            ft_error('Some electrodes are too distant from the plane: consider repositioning them')
          else
            % project point on plane
            Ppr = pointproj(P,[headmodel.pos headmodel.ori]);
            pos(i,:) = Ppr;
          end
        end
      end
      sens.elecpos = pos;
      
    case {'slab_monopole'}
      % electrodes' all-to-all distances
      numelc  = size(sens.elecpos,1);
      ref_elc = sens.elecpos(1,:);
      md  = dist( (sens.elecpos-repmat(ref_elc,[numelc 1]))' );
      % choose min distance between electrodes
      md  = min(md(1,2:end));
      pos = sens.elecpos;
      % looks for contacts outside the strip which are not too far away
      % and projects them on the nearest plane
      for i=1:size(pos,1)
        P = pos(i,:);
        instrip1 = acos(dot(headmodel.ori1,(P-headmodel.pos1)./norm(P-headmodel.pos1))) > pi/2;
        instrip2 = acos(dot(headmodel.ori2,(P-headmodel.pos2)./norm(P-headmodel.pos2))) > pi/2;
        is_in_empty = ~(instrip1&instrip2);
        if is_in_empty
          dPplane1 = abs(dot(headmodel.ori1, headmodel.pos1-P, 2));
          dPplane2 = abs(dot(headmodel.ori2, headmodel.pos2-P, 2));
          if dPplane1>md && dPplane2>md
            ft_error('Some electrodes are too distant from the planes: consider repositioning them')
          elseif dPplane2>dPplane1
            % project point on nearest plane
            Ppr = pointproj(P,[headmodel.pos1 headmodel.ori1]);
            pos(i,:) = Ppr;
          else
            % project point on nearest plane
            Ppr = pointproj(P,[headmodel.pos2 headmodel.ori2]);
            pos(i,:) = Ppr;
          end
        end
      end
      sens.elecpos = pos;
      
    case {'singlesphere', 'concentricspheres'}
      % ensure that the electrodes ly on the skin surface
      radius = max(headmodel.r);
      pos    = sens.elecpos;
      if isfield(headmodel, 'o')
        % shift the the centre of the sphere to the origin
        pos(:,1) = pos(:,1) - headmodel.o(1);
        pos(:,2) = pos(:,2) - headmodel.o(2);
        pos(:,3) = pos(:,3) - headmodel.o(3);
      end
      distance = sqrt(sum(pos.^2,2)); % to the center of the sphere
      if any((abs(distance-radius)/radius)>0.005)
        ft_warning('electrodes do not lie on skin surface -> using radial projection')
      end
      pos = pos * radius ./ [distance distance distance];
      if isfield(headmodel, 'o')
        % shift the center back to the original location
        pos(:,1) = pos(:,1) + headmodel.o(1);
        pos(:,2) = pos(:,2) + headmodel.o(2);
        pos(:,3) = pos(:,3) + headmodel.o(3);
      end
      sens.elecpos = pos;
      
    case {'bem', 'dipoli', 'asa', 'bemcp'}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % do postprocessing of volume and electrodes in case of BEM model
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % project the electrodes on the skin and determine the bilinear interpolation matrix
      if ~isfield(headmodel, 'tra') && (isfield(headmodel, 'mat') && ~isempty(headmodel.mat))
        % determine boundary corresponding with skin and inner_skull_surface
        if ~isfield(headmodel, 'skin_surface')
          headmodel.skin_surface = find_outermost_boundary(headmodel.bnd);
          fprintf('determining skin compartment (%d)\n', headmodel.skin_surface);
        end
        if ~isfield(headmodel, 'source')
          headmodel.source = find_innermost_boundary(headmodel.bnd);
          fprintf('determining source compartment (%d)\n', headmodel.source);
        end
        if size(headmodel.mat,1)~=size(headmodel.mat,2) && size(headmodel.mat,1)==length(sens.elecpos)
          fprintf('electrode transfer and system matrix were already combined\n');
        else
          fprintf('projecting electrodes on skin surface\n');
          % compute linear interpolation from triangle vertices towards electrodes
          [el, prj] = project_elec(sens.elecpos, headmodel.bnd(headmodel.skin_surface).pos, headmodel.bnd(headmodel.skin_surface).tri);
          tra       = transfer_elec(headmodel.bnd(headmodel.skin_surface).pos, headmodel.bnd(headmodel.skin_surface).tri, el);
          
          % replace the original electrode positions by the projected positions
          sens.elecpos = prj;
          
          if size(headmodel.mat,1)==size(headmodel.bnd(headmodel.skin_surface).pos,1)
            % construct the transfer from only the skin vertices towards electrodes
            interp = tra;
          else
            % construct the transfer from all vertices (also inner_skull_surface/outer_skull_surface) towards electrodes
            interp = [];
            for i=1:length(headmodel.bnd)
              if i==headmodel.skin_surface
                interp = [interp, tra];
              else
                interp = [interp, zeros(size(el,1), size(headmodel.bnd(i).pos,1))];
              end
            end
          end
          
          % incorporate the linear interpolation matrix and the system matrix into one matrix
          % this speeds up the subsequent repeated leadfield computations
          fprintf('combining electrode transfer and system matrix\n');
          
          % convert to sparse matrix to speed up the subsequent multiplication
          interp  = sparse(interp);
          headmodel.mat = interp * headmodel.mat;
          % ensure that the model potential will be average referenced
          avg = mean(headmodel.mat, 1);
          headmodel.mat = headmodel.mat - repmat(avg, size(headmodel.mat,1), 1);
        end
      end
    case  'openmeeg' 
        % don't do anything, h2em or h2mm generated later in ft_prepare_leadfield
      
    case 'fns'
      if isfield(headmodel,'bnd')
        [el, prj] = project_elec(sens.elecpos, headmodel.bnd.pos, headmodel.bnd.tri);
        sens.tra = transfer_elec(headmodel.bnd.pos, headmodel.bnd.tri, el);
        % replace the original electrode positions by the projected positions
        sens.elecpos = prj;
      end
      
    case 'simbio'
      % check that the external toolbox is present
      ft_hastoolbox('simbio', 1);
      
      % extract the outer surface
      bnd = mesh2edge(headmodel);
      for j=1:length(sens.label)
        d = bsxfun(@minus, bnd.pos, sens.elecpos(j,:));
        [d, i] = min(sum(d.^2, 2));
        % replace the position of each electrode by the closest vertex
        sens.elecpos(j,:) = bnd.pos(i,:);
      end
      
      headmodel.transfer = sb_transfer(headmodel,sens);
      
    case 'interpolate'
      % this is to allow moving leadfield files
      if ~exist(headmodel.filename{1}, 'file')
        for i = 1:length(headmodel.filename)
          [p, f, x] = fileparts(headmodel.filename{i});
          headmodel.filename{i} = fullfile(vpath, [f x]);
        end
      end
      
      matchlab = isequal(sens.label, headmodel.sens.label);
      matchpos = isequal(sens.elecpos, headmodel.sens.elecpos);
      matchtra = (~isfield(sens, 'tra') && ~isfield(headmodel.sens, 'tra')) || isequal(sens.tra, headmodel.sens.tra); 

      if matchlab && matchpos && matchtra
        % the input sensor array matches precisely with the forward model
        % no further interpolation is needed
      else
        % interpolate the channels in the forward model to the desired channels
        filename = tempname;
        headmodel  = ft_headmodel_interpolate(filename, sens, headmodel);
        % update the sensor array with the one from the volume conductor
        sens = headmodel.sens;
      end % if recomputing interpolation
      
      % for the leadfield computations the @nifti object is used to map the image data into memory
      ft_hastoolbox('spm8up', 1);
      for i=1:length(headmodel.sens.label)
        % map each of the leadfield files into memory
        headmodel.chan{i} = nifti(headmodel.filename{i});
      end
      
    otherwise
      ft_error('unsupported volume conductor model for EEG');
  end
  
  % FIXME this needs careful thought to ensure that the average referencing which is now done here and there, and that the linear interpolation in case of BEM are all dealt with consistently
  % % always ensure that there is a linear transfer matrix for
  % % rereferencing the EEG potential
  % if ~isfield(sens, 'tra');
  %   sens.tra = eye(length(sens.label));
  % end
  
  % update the channel positions as the electrodes were projected to the skin surface
  [pos, ori, lab] = channelposition(sens);
  [selsens, selpos] = match_str(sens.label, lab);
  sens.chanpos = nan(length(sens.label),3);
  sens.chanpos(selsens,:) = pos(selpos,:);

end % if iseeg or ismeg

if isfield(sens, 'tra')
  if issparse(sens.tra) && size(sens.tra, 1)==1
    % this multiplication would result in a sparse leadfield, which is not what we want
    % the effect can be demonstrated as sparse(1)*rand(1,10), see also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1169#c7
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
