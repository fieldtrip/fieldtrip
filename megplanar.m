function [interp] = megplanar(cfg, data);

% MEGPLANAR computes planar MEG gradients gradients for raw data
% obtained from PREPROCESSING or an average ERF that was computed using
% TIMELOCKANALYSIS. It can also work on data in the frequency domain, 
% obtained with FREQANALYSIS. Prerequisite for this is that the data contain
% complex-valued fourierspectra.
%
% Use as
%    [interp] = megplanar(cfg, data)
%
% The configuration should contain
%   cfg.planarmethod   = 'orig' | 'sincos' | 'fitplane' | 'sourceproject'
%   cfg.channel        =  Nx1 cell-array with selection of channels (default = 'MEG'),
%                         see CHANNELSELECTION for details
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%
% The methods orig, sincos and fitplane are all based on a neighbourhood
% interpolation. For these methods you can specify
%   cfg.neighbourdist  = default is 4 cm 
% 
% In the 'sourceproject' method a minumum current estimate is done using a
% large number of dipoles that are placed in the upper layer of the brain
% surface, followed by a forward computation towards a planar gradiometer
% array. This requires the specification of a volume conduction model of
% the head and of a source model. The 'sourceproject' method is not supported for
% frequency domain data.
%
% A head model must be specified with
%   cfg.hdmfile     = string, file containing the volume conduction model
% or alternatively manually using 
%   cfg.vol.r       = radius of sphere
%   cfg.vol.o       = [x, y, z] position of origin
%
% A dipole layer representing the brain surface must be specified with
%   cfg.inwardshift = depth of the source layer relative to the head model surface (default = 2.5, which is adequate for a skin-based head model)
%   cfg.spheremesh  = number of dipoles in the source layer (default = 642)
%   cfg.pruneratio  = for singular values, default is 1e-3
%   cfg.headshape   = a filename containing headshape, a structure containing a
%                     single triangulated boundary, or a Nx3 matrix with surface
%                     points
% If no headshape is specified, the dipole layer will be based on the inner compartment
% of the volume conduction model.
% 
% See also COMBINEPLANAR

% This function depends on PREPARE_BRAIN_SURFACE which has the following options:
% cfg.headshape  (default set in MEGPLANAR: cfg.headshape = 'headmodel'), documented
% cfg.inwardshift (default set in MEGPLANAR: cfg.inwardshift = 2.5), documented
% cfg.spheremesh (default set in MEGPLANAR: cfg.spheremesh = 642), documented
%
% This function depends on PREPARE_VOL_SENS which has the following options:
% cfg.channel
% cfg.elec
% cfg.elecfile
% cfg.grad
% cfg.gradfile
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: megplanar.m,v $
% Revision 1.45  2009/07/19 13:19:42  jansch
% fixed re-conversion into timelock datatype when input is timelock datatype
%
% Revision 1.44  2009/07/17 08:05:12  jansch
% fixed small bug
%
% Revision 1.43  2009/07/16 15:07:26  jansch
% updated documentation
%
% Revision 1.42  2009/07/16 15:05:15  jansch
% added possibility to perform planar transformation on frequency domain data.
% the frequency data of course should contain fourierspectra
%
% Revision 1.41  2009/07/14 12:35:36  roboos
% ensure that the input data is meg with a grad field
%
% Revision 1.40  2009/06/04 10:03:46  roboos
% fixed problem in the input og cfg.headshape when it was not required
%
% Revision 1.39  2009/05/25 08:05:18  roboos
% ensure that cfg.headshape is a sturct and not a config object (in case tracking is on)
%
% Revision 1.38  2009/05/14 19:23:03  roboos
% consistent handling of cfg.headshape in code and documentation
%
% Revision 1.37  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.36  2008/11/21 12:48:17  sashae
% added call to checkconfig at start and end of function
%
% Revision 1.35  2008/11/13 21:53:57  roboos
% ensure that grad.unit is specified, also support bti148
%
% Revision 1.34  2008/10/02 15:32:20  sashae
% replaced call to createsubcfg with checkconfig
%
% Revision 1.33  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.32  2008/07/15 19:56:44  roboos
% moved cfg details for dipole grid to subcfg (cfg.grid)subcfg (cfg.grid.xxx)
%
% Revision 1.31  2008/05/08 11:12:19  jansch
% added support for bti248 magnetometer system
%
% Revision 1.30  2008/05/06 15:43:46  sashae
% change in trial selection, cfg.trials can be logical
%
% Revision 1.29  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.28  2008/01/31 17:20:16  sashae
% added option for trial selection
%
% Revision 1.27  2007/05/30 13:22:43  roboos
% added cfg.channel default to ensure that only MEG channels (and not reference channels) are converted to planar
%
% Revision 1.26  2007/05/29 14:26:44  ingnie
% changed sensortype to senstype (in calling checkdata) to avoid overlap with function name sensortype
%
% Revision 1.25  2007/05/29 12:51:31  roboos
% added new options for checkdata
%
% Revision 1.24  2007/05/02 15:59:13  roboos
% be more strict on the input and output data: It is now the task of
% the private/checkdata function to convert the input data to raw
% data (i.e. as if it were coming straight from preprocessing).
% Furthermore, the output data is NOT converted back any more to the
% input data, i.e. the output data is the same as what it would be
% on raw data as input, regardless of the actual input.
%
% Revision 1.23  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.22  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.21  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.20  2006/11/24 11:14:08  roboos
% changed default planarmethod from orig into sincos
%
% Revision 1.19  2006/07/24 07:59:16  roboos
% updated documentation
%
% Revision 1.18  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.17  2006/04/06 16:17:10  ingnie
% updated documentation
%
% Revision 1.16  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.15  2006/01/30 14:34:59  roboos
% added two comments as a reminder that prepare_dipole_grid should be used in the future
%
% Revision 1.14  2006/01/30 14:22:11  roboos
% renamed axial2planar into constructplanargrad
%
% Revision 1.13  2005/12/14 14:16:07  roboos
% switched from prepare_brain_surface to headsurface subfunction
%
% Revision 1.12  2005/07/29 07:18:52  roboos
% ensure that non-gradiometer channels are also handled correctly for planarmethod=sourceproject
%
% Revision 1.11  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.10  2005/06/08 13:40:33  roboos
% replaced the specific call to either meg_leadfield or eeg_leadfield to the generic compute_leadfield
%
% Revision 1.9  2005/06/02 12:18:41  roboos
% changed handling of input data: All input data that contains averages is converted to raw trials (like the output from preprocessing) prior to further processing. The output data is converted back into a format similar to the original input data using RAW2MEG.
%
% Revision 1.8  2004/05/07 14:13:50  roberto
% small change in documentation
%
% Revision 1.7  2004/05/05 13:55:20  roberto
% fixed bug that occured when average data was made with keeptrials=yes
%
% Revision 1.6  2004/04/28 10:05:49  roberto
% added planar gradiometer definition to output (in case of nearest neighbours approach)
%
% Revision 1.5  2004/04/13 16:31:09  roberto
% fixed bug in dbstack selection of function filename for Matlab 6.1
%
% Revision 1.4  2004/04/13 14:25:24  roberto
% wrapped code-snippet around mfilename to make it compatible with Matlab 6.1
%
% Revision 1.3  2004/04/08 15:49:21  roberto
% in the sourceproject method removed the selection of channels, now using prepare_vol_sens.
% replaced the construction of the brain dipole surface with a call to the new function
% prepare_brain_surface (shared with meg_realign)
%
% Revision 1.2  2004/04/06 20:14:35  roberto
% Note: the previous revision was not merely a copy-and-paste version
% of meginterpolate. It also added support for different planar
% gradient computation methods, including the original by Ole and it
% extended all methods to be able to work with less than 151 channels.
% This revision has only changes in help and comments.
%
% Revision 1.1  2004/04/06 20:01:49  roberto
% created new implementation by moving all functionality from meginterpolate to this function
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

isfreq = datatype(data, 'freq');
israw  = datatype(data, 'raw');
istlck = datatype(data, 'timelock');

% check if the input data is valid for this function
data  = checkdata(data, 'datatype', {'raw' 'freq'}, 'feedback', 'yes', 'ismeg', 'yes', 'senstype', {'ctf151', 'ctf275', 'bti148', 'bti248'});

if isfreq,
  if ~isfield(data, 'fourierspctrm'), error('freq data should containt fourier spectra'); end
end

% set the default configuration 
if ~isfield(cfg, 'channel'),       cfg.channel = 'MEG';             end
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';              end
if ~isfield(cfg, 'neighbourdist'), cfg.neighbourdist = 4;           end
if ~isfield(cfg, 'planarmethod'),  cfg.planarmethod = 'sincos';     end
if strcmp(cfg.planarmethod, 'sourceproject')
  if ~isfield(cfg, 'headshape'),     cfg.headshape = [];            end % empty will result in the vol being used
  if ~isfield(cfg, 'inwardshift'),   cfg.inwardshift = 2.5;         end
  if ~isfield(cfg, 'pruneratio'),    cfg.pruneratio = 1e-3;         end
  if ~isfield(cfg, 'spheremesh'),    cfg.spheremesh = 642;          end
end

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% put the low-level options pertaining to the dipole grid in their own field
cfg = checkconfig(cfg, 'createsubcfg',  {'grid'});
cfg = checkconfig(cfg, 'renamedvalue',  {'headshape', 'headmodel', []});

% select trials of interest
if ~strcmp(cfg.trials, 'all') && israw
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
  % update the trial definition (trl)
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    trl = findcfg(data.cfg, 'trl');
  else
    trl = [];
  end
  if isempty(trl)
    % a trial definition is expected in each continuous data set
    warning('could not locate the trial definition ''trl'' in the data structure');
  else
    cfg.trlold=trl;
    cfg.trl=trl(cfg.trials,:);
  end
elseif ~strcmp(cfg.trials, 'all') && isfreq
  warning('subselection of trials is only supported for raw data as input');
end

if israw || istlck,
  Ntrials = length(data.trial);
elseif isfreq,
  Ntrials = length(data.cumtapcnt);
  Nfreq   = length(data.freq);
  if isfield(data, 'time')
    Ntime = length(data.time);
  else
    Ntime = 1;
  end
end
Nchan   = length(data.label);

% find the corresponding channels in the data and the gradiometer array
cfg.channel = channelselection(cfg.channel, data.label);
selindx = match_str(data.label, cfg.channel);  % these are selected to be transformed into planar
[gradindx, dataindx] = match_str(data.grad.label, data.label(selindx));
dataindx = selindx(dataindx); % re-index the data channel selection
fprintf('computing planar gradient for %d channels\n', length(dataindx));

% ensure that they are row-vectors
gradindx = gradindx(:)';
dataindx = dataindx(:)';

% make a copy for convenience
grad = data.grad;
% ensure that the units are specified
grad = convert_units(grad);
% remove the gradiometers that have no channel in the data
% we only need the bottom coil for the position of the channel
grad.label = grad.label(gradindx);
grad.pnt   = grad.pnt(gradindx, :);
grad.ori   = grad.ori(gradindx, :);
Ngrad = length(gradindx);

for i=1:Ngrad
  for j=(i+1):Ngrad  
    distance(i,j) = norm(grad.pnt(i,:)-grad.pnt(j,:));
    distance(j,i) = distance(i,j);
  end
end

fprintf('minimum distance between gradiometers is %6.2f %s\n', min(distance(find(distance~=0))), grad.unit);
fprintf('maximum distance between gradiometers is %6.2f %s\n', max(distance(find(distance~=0))), grad.unit);

% select the channels that are neighbours, channel is not a neighbour of itself
neighbsel = distance<cfg.neighbourdist;
for i=1:Ngrad
  neighbsel(i,i) = 0;
end
fprintf('average number of neighbours is %f\n', sum(neighbsel(:))./size(neighbsel,1));

gradH = zeros(Ngrad, Ngrad);
gradV = zeros(Ngrad, Ngrad);

if strcmp(cfg.planarmethod, 'orig')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This is the original method from Ole.  It has a different way of
  % making the coordinate transformation that I do not fully understand.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % get the positions of bottom and top coil
  X  = grad.pnt(:, 1);
  Y  = grad.pnt(:, 2);
  Z  = grad.pnt(:, 3);
  X2 = grad.pnt(:, 1) + grad.ori(:,1);
  Y2 = grad.pnt(:, 2) + grad.ori(:,2);
  Z2 = grad.pnt(:, 3) + grad.ori(:,3);
  
  for k=1:Ngrad
    % translate the current coil to the origin
    Xc = X - X(k);
    Yc = Y - Y(k);
    Zc = Z - Z(k);
    X2c = X2 - X(k);
    Y2c = Y2 - Y(k);
    Z2c = Z2 - Z(k);
    
    X = Xc;
    Y = Yc;
    Z = Zc;
    X2 = X2c;
    Y2 = Y2c;
    Z2 = Z2c;
    
    % rotate around z-axis
    PhiZ = -atan(Y2(k)/(0.00000001+X2(k)));
    Xc  = X*cos(PhiZ) - Y*sin(PhiZ);
    Yc  = X*sin(PhiZ) + Y*cos(PhiZ);
    Zc  = Z;
    X2c = X2*cos(PhiZ) - Y2*sin(PhiZ);
    Y2c = X2*sin(PhiZ) + Y2*cos(PhiZ);
    Z2c = Z2;
    
    X = Xc;
    Y = Yc;
    Z = Zc;
    X2 = X2c;
    Y2 = Y2c;
    Z2 = Z2c;
    
    % rotate around y-axis
    PhiY = atan(Z2(k)/(0.00000001+X2(k)));
    Xc  = X*cos(PhiY) + Z*sin(PhiY);
    Yc  = Y;
    Zc  = -X*sin(PhiY) + Z*cos(PhiY);
    X2c = X2*cos(PhiY) + Z2*sin(PhiY);
    Y2c = Y2;
    Z2c = -X2*sin(PhiY) + Z2*cos(PhiY);;
    
    X = Xc;
    Y = Yc;
    Z = Zc;
    X2 = X2c;
    Y2 = Y2c;
    Z2 = Z2c;
    
    % compute planar gradients in y and z directions
    for l=find(neighbsel(k,:))
        ac = Y(l)/distance(k,l);
        gradH(l,k) = ac/distance(k,l);
        gradH(l,l) = gradH(l,l) - ac/distance(k,l);
        
        ac = Z(l)/distance(k,l);
        gradV(l,k) = ac/distance(k,l);
        gradV(l,l) = gradV(l,l) - ac/distance(k,l);    
    end % for l
  end % for k
  
elseif strcmp(cfg.planarmethod, 'sincos')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This attempts to re-implements Ole's method, exept that the definition of the 
  % horizontal and vertical direction is different.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for chan=1:Ngrad
    % Attach a local coordinate system to this gradiometer:
    % the origin at the location of its bottom coil
    % the z-axis pointing outwards from the head
    % the x-axis pointing horizontal w.r.t. the head
    % the y-axis pointing vertical, i.e. approximately towards the vertex
    this_o = grad.pnt(chan,:);
    this_z = grad.ori(chan,:);          this_z = this_z / norm(this_z);
    this_x = cross([0 0 1], this_z);
    if all(this_x==0)
      this_x = [1 0 0];
    else
      this_x = this_x / norm(this_x);
    end
    this_y = cross(this_z, this_x);
    
    for neighb=find(neighbsel(chan, :))
      vec = grad.pnt(neighb,:) - this_o;    % vector from sensor to neighbour
      proj_x = dot(vec, this_x);            % projection along x-axis (horizontal)
      proj_y = dot(vec, this_y);            % projection along y-axiz (vertical)
      proj_z = dot(vec, this_z);            % projection along z-axis 
      
      gradH(chan, chan)   = gradH(chan,chan)    - proj_x / (norm(vec).^2);
      gradH(chan, neighb) =                       proj_x / (norm(vec).^2);
      gradV(chan, chan)   = gradV(chan,chan)    - proj_y / (norm(vec).^2);
      gradV(chan, neighb) =                       proj_y / (norm(vec).^2);
    end
  end
  
elseif strcmp(cfg.planarmethod, 'fitplane')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fit a plane through the B=f(x,y) plane and compute its two gradients
  % The first point in the plane is the gradiometer itself,
  % the neighbours are the subsequent points. This method also returns the
  % offset of the B-plane at each sensor, which is appriximately equal to the
  % field itself.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for chan=1:Ngrad
    % Attach a local coordinate system to this gradiometer:
    % the origin at the location of its bottom coil
    % the z-axis pointing outwards from the head
    % the x-axis pointing horizontal w.r.t. the head
    % the y-axis pointing vertical, i.e. approximately towards the vertex
    this_o = grad.pnt(chan,:);
    this_z = grad.ori(chan,:);         
    this_z = this_z / norm(this_z);
    this_x = cross([0 0 1], this_z);
    if all(this_x==0)
      this_x = [1 0 0];
    else
      this_x = this_x / norm(this_x);
    end
    this_y = cross(this_z, this_x);
    
    % add the relative position in local coordinates to the list of points
    % starting with the channel position itself (which is the local origin)
    x = 0;
    y = 0;
    e = 1;
    neighbindx = find(neighbsel(chan, :));
    for neighb=neighbindx
      vec = grad.pnt(neighb,:) - this_o;          % vector from sensor to neighbour
      x(end+1) = dot(vec, this_x);                % projection along x-axis (horizontal)
      y(end+1) = dot(vec, this_y);                % projection along y-axiz (vertical)
      e(end+1) = 1;                               % required to fit the constant
    end
    A = pinv([x(:) y(:) e(:)]);
    gradH(chan,[chan neighbindx]) = A(1,:);
    gradV(chan,[chan neighbindx]) = A(2,:);
    gradC(chan,[chan neighbindx]) = A(3,:);
  end
  
elseif strcmp(cfg.planarmethod, 'sourceproject') && israw
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Do an inverse computation with a simplified distributed source model 
  % and compute forward again with the axial gradiometer array replaced by
  % a planar one.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % PREPARE_VOL_SENS will match the data labels, the gradiometer labels and the 
  % volume model labels (in case of a multisphere model) and result in a gradiometer 
  % definition that only contains the gradiometers that are present in the data.
  [vol, axial.grad, cfg] = prepare_headmodel(cfg, data);

  % determine the dipole layer that represents the surface of the brain
  if isempty(cfg.headshape)
    % construct from the inner layer of the volume conduction model
    pos = headsurface(vol, axial.grad, 'surface', 'cortex', 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  else
    % get the surface describing the head shape
    if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
      % use the headshape surface specified in the configuration
      headshape = cfg.headshape;
    elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
      % use the headshape points specified in the configuration
      headshape.pnt = cfg.headshape;
    elseif ischar(cfg.headshape)
      % read the headshape from file
      headshape = read_headshape(cfg.headshape);
    else
      error('cfg.headshape is not specified correctly')
    end
    if ~isfield(headshape, 'tri')
      % generate a closed triangulation from the surface points
      headshape.pnt = unique(headshape.pnt, 'rows');
      headshape.tri = projecttri(headshape.pnt);
    end
    % construct from the head surface
    pos = headsurface([], [], 'headshape', headshape, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  end

  % compute the forward model for the axial gradiometers
  fprintf('computing forward model for %d dipoles\n', size(pos,1));
  lfold = compute_leadfield(pos, axial.grad, vol);

  % construct the planar gradient definition and compute its forward model
  % this will not work for a multisphere model, compute_leadfield will catch
  % the error
  planar.grad = constructplanargrad([], axial.grad);
  lfnew = compute_leadfield(pos, planar.grad, vol);

  % compute the interpolation matrix
  transform = lfnew * prunedinv(lfold, cfg.pruneratio);

  % interpolate the data towards the planar gradiometers
  for i=1:Ntrials
    fprintf('interpolating trial %d to planar gradiometer\n', i);
    interp.trial{i} = transform * data.trial{i}(dataindx,:);
  end % for Ntrials

  % all planar gradiometer channels are included in the output
  interp.grad  = planar.grad;
  interp.label = planar.grad.label;

  % copy the non-gradiometer channels back into the output data
  other = setdiff(1:Nchan, dataindx);
  for i=other
    interp.label{end+1} = data.label{i};
    for j=1:Ntrials
      interp.trial{j}(end+1,:) = data.trial{j}(i,:);
    end
  end

elseif strcmp(cfg.planarmethod, 'sourceproject') && isfreq
  error('the method ''sourceproject'' is not supported for frequency data as input');
else
  error('unknown method for computation of planar gradient');
end % cfg.planarmethod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the linear transformation for any of the neighbourhood methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(cfg.planarmethod, 'sourceproject')
  % these exist always, except in the sourceproject method
  transformH = zeros(Ngrad,Nchan);
  transformV = zeros(Ngrad,Nchan);
  % each row of the transformation matrix should correspond with the order
  % of the channels in the original data
  transformH(:,dataindx) = gradH;
  transformV(:,dataindx) = gradV;
  
  if strcmp(cfg.planarmethod, 'fitplane')
    % this only exists if the fitplane method was used
    transformC = zeros(Ngrad,Nchan);
    transformC(:,dataindx) = gradC;
  else
    transformC = [];
  end
  
  other = setdiff(1:Nchan, dataindx);          % other channels not having a planar gradient, e.g. EMG
  transformO = zeros(length(other), Nchan);
  for chan=1:length(other)
    transformO(chan,other(chan)) = 1;          % these should stay in the data with a weight of one
  end
  
  % rename the labels to match the new channel content
  labelH = {};
  labelV = {};
  labelC = {};
  labelO = {};
  for i=dataindx
    labelH{end+1} = sprintf('%s_dH', data.label{i});
  end
  for i=dataindx
    labelV{end+1} = sprintf('%s_dV', data.label{i});
  end
  if strcmp(cfg.planarmethod, 'fitplane')
    for i=dataindx
      labelC{end+1} = sprintf('%s_C', data.label{i});
    end
  end
  for i=other
    % these keep the original label
    labelO{end+1} = data.label{i};
  end
  
  % combine the different sections in one linear transformation matrix
  transform = [transformH; transformV; transformC; transformO];

  % combine the new labels into a single cell-array
  interp.label = [labelH(:); labelV(:); labelC(:); labelO(:)];

  % construct the planar gradiometer definition
  planar.grad.pnt = data.grad.pnt;		% coils are on the same location
  planar.grad.ori = data.grad.ori;		% coils have the same orientation
  planar.grad.tra = [gradH; gradV] * data.grad.tra(gradindx,:); 
  planar.grad.label = [labelH(:); labelV(:)];
  try,
    planar.grad.unit = data.grad.unit;
  end

  % remember the planar gradiometer definition
  interp.grad = planar.grad;
  
  if israw || istlck,
    % compute the planar gradient by multiplying the gradient-matrices with the data
    for trial=1:Ntrials
      interp.trial{trial} = transform * data.trial{trial};
    end
  
    % these should be remembered from the original data
    interp.fsample = data.fsample;
    interp.time    = data.time;
  
    if istlck
      % convert back into timelock structure
      interp = checkdata(interp, 'datatype', 'timelock');
    end

  elseif isfreq
    % compute the planar gradient by multiplying the gradient-matrices with the data
    siz       = size(data.fourierspctrm);
    planardat = zeros(siz(1), size(transform,1), Nfreq, Ntime);
    for foilop=1:Nfreq
      for timlop = 1:Ntime
        planardat(:,:,foilop,timlop) = data.fourierspctrm(:,:,foilop,timlop) * transform';
      end
    end
    interp.fourierspctrm = planardat;

    if isfield(data, 'time'), interp.time = data.time; end
    if isfield(data, 'cumtapcnt'), interp.cumtapcnt = data.cumtapcnt; end
    if isfield(data, 'cumsumcnt'), interp.cumsumcnt = data.cumsumcnt; end
    interp.freq   = data.freq;
    interp.dimord = data.dimord;  
    interp.transform = transform;
  end

end % nearest neighbours methods

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% store the configuration of this function call, including that of the previous function call
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: megplanar.m,v 1.45 2009/07/19 13:19:42 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
interp.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the inverse using a pruned SVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lfi] = prunedinv(lf, r)
[u, s, v] = svd(lf);
p = find(s<(s(1,1)*r) & s~=0);
fprintf('pruning %d out of %d singular values\n', length(p), min(size(s)));
s(p) = 0;
s(find(s~=0)) = 1./s(find(s~=0));
lfi = v * s' * u';
