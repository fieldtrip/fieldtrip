function vol = ft_headmodel_localspheres(geometry, grad, varargin)
  
% FT_HEADMODEL_LOCALSPHERES constructs a MEG volume conduction model in
% with a local sphere fitted to the head or brain surface for each separate
% channel
%
% This implements
%   Huang MX, Mosher JC, Leahy RM. "A sensor-weighted overlapping-sphere
%   head model and exhaustive head model comparison for MEG." Phys Med
%   Biol. 1999 Feb;44(2):423-40
%
% Use as
%   vol = ft_headmodel_localspheres(geom, grad, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   feedback = boolean, true or false
%   radius   = number, radius of sphere within which headshape points will
%     be included for the fitting algorithm
%   maxradius = number, if for a given sensor the fitted radius exceeds
%     this value, the radius and origin will be replaced witht the single
%     sphere fit
%   baseline  = number
%
% See also FT_PREPARE_HEADMODEL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% get the additional inputs and set the defaults
% headshape    = keyval('headshape', varargin);
feedback     = keyval('feedback',  varargin); if isempty(feedback), feedback = true; end
radius       = keyval('radius',    varargin); if isempty(radius),   radius   = 8.5;  end
maxradius    = keyval('maxradius', varargin); if isempty(maxradius),maxradius = 20;  end
baseline     = keyval('baseline',  varargin); if isempty(baseline), baseline  = 5;   end

if ischar(feedback)
  if strcmp(feedback, 'yes') 
    feedback = true;
  else 
    feedback = false;
  end
end

% start with an empty volume conductor
vol = [];

if isfield(geometry, 'unit')
  % copy the geometrical units into he volume conductor
  vol.unit = geometry.unit;
end

Nshape = size(geometry.pnt,1);
Nchan  = numel(grad.label);

% set up an empty figure
if istrue(feedback)
  clf
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
end

% plot all channels and headshape points
if istrue(feedback)
  cla
  ft_plot_sens(grad);
  ft_plot_mesh(geometry, 'vertexcolor', 'g', 'facecolor', 'none', 'edgecolor', 'none');
  drawnow
end

% fit a single sphere to all headshape points
[single_o, single_r] = fitsphere(geometry.pnt);
fprintf('single sphere,   %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', Nshape, single_o(1), single_o(2), single_o(3), single_r);


% allocate empty matrices that will hold the results
vol.r = zeros(Nchan,1);    % radius of every sphere
vol.o = zeros(Nchan,3);    % origin of every sphere
vol.label = cell(Nchan,1); % corresponding gradiometer channel label for every sphere

for chan=1:Nchan
  coilsel = find(grad.tra(chan,:)~=0);
  allpnt  = grad.pnt(coilsel, :);   % position of all coils belonging to this channel
  allori  = grad.ori(coilsel, :);   % orientation of all coils belonging to this channel
  
  if istrue(feedback)
    cla
    plot3(grad.pnt(:,1), grad.pnt(:,2), grad.pnt(:,3), 'b.');   % all coils
    plot3(allpnt(:,1), allpnt(:,2), allpnt(:,3), 'r*');     % this channel in red
  end
  
  % determine the average position and orientation of this channel
  thispnt = mean(allpnt,1);
  [u, s, v] = svd(allori);
  thisori = v(:,1)';
  if dot(thispnt,thisori)<0
    % the orientation should be outwards pointing
    thisori = -thisori;
  end
  
  % compute the distance from every coil along this channels orientation
  dist = zeros(size(coilsel));
  for i=1:length(coilsel)
    dist(i) = dot((allpnt(i,:)-thispnt), thisori);
  end
  
  [m, i] = min(dist);
  % check whether the minimum difference is larger than a typical distance
  if abs(m)>(baseline/4)
    % replace the position of this channel by the coil that is the closest to the head (axial gradiometer)
    % except when the center of the channel is approximately just as good (planar gradiometer)
    thispnt = allpnt(i,:);
  end
  
  % find the headshape points that are close to this channel
  dist = sqrt(sum((geometry.pnt-repmat(thispnt,Nshape,1)).^2, 2));
  shapesel = find(dist<radius);
  if feedback
    ft_plot_mesh(geometry.pnt(shapesel,:), 'vertexcolor', 'g');
    drawnow
  end
  
  % fit a sphere to these headshape points
  if length(shapesel)>10
    [o, r] = fitsphere(geometry.pnt(shapesel,:));
    fprintf('channel = %s, %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', grad.label{chan}, length(shapesel), o(1), o(2), o(3), r);
  else
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end
  
  if r > maxradius
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end
  
  % add this sphere to the volume conductor
  vol.o(chan,:)   = o;
  vol.r(chan)     = r;
  vol.label{chan} = grad.label{chan};
end % for all channels

vol.type = 'multisphere';

