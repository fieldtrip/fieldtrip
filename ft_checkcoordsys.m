function [varargout] = ft_checkcoordsys(cfg, varargin)

% FT_CHECKCOORDSYS allows the user t o perform  a check on the coordinatesystem
% and units in which the input object is expressed. 
%
% Use as
%   [output] = ft_checkcoordsys(cfg, data)
% or as
%   [output1, output2, ...] = ft_checkcoordsys(cfg, data1, data2, ...)
%
% 

ndata   = numel(varargin);
nargout = ndata;
for k = 1:ndata
  dtype{k} = datatype(varargin{k});

  % try to ensure the data to have a unit FIXME consider to move this to checkdata
  if k==1,
    try,
      varargin{k} = ft_convert_units(varargin{k});
    end
    unit = varargin{k}.unit;
  else
    try,
      varargin{k} = ft_convert_units(varargin{k}, unit);
    end
  end
end

if exist('unit', 'var')
  switch unit
  case 'mm'
    axmax = 150;
    rbol  = 5;
  case 'cm'
    axmax = 15;
    rbol  = 0.5;
  case 'm'
    axmax = 0.15;
    rbol  = 0.005;
  otherwise
    % assume cm
    axmax = 15;
    rbol  = 5;
  end
else
  axmax = 15;
end

for k = 1:ndata
  
  % do some checks on the geometrical object
  data = varargin{k};
  switch dtype{k}
  case 'volume'
    % contains transform
    transform    = data.transform;
    hastransform = 1;
    haspos       = 0;
  case 'source'
    % contains pos
    pos          = data.pos;
    hastransform = 0;
    haspos       = 1;
    hastri       = 0; %FIXME may need to be supported in the future
  case 'dip'
    % contains pos
    pos          = data.pos;
    hastransform = 0;
    haspos       = 1;
    hastri       = 0; %FIXME may need to be supported in the future
  case {'freq', 'raw', 'timelock', 'spike', 'mvar', 'freqmvar'}
    % these may contain geometrical objects in grad or elec fields
    if isfield(data, 'grad')
      tmpdata = data.grad;
      tmpdata = ft_checkcoordinates([], tmpdata);
      data.grad = tmpdata;
    elseif isfield(data, 'elec')
      tmpdata = data.elec;
      tmpdata = ft_checkcoordinates([], tmpdata);
      data.grad = tmpdata;
    end
    %haspos       = 1;
  case 'unknown'
    % other geometrical objects typically contain geo info in pnt-field
    if isfield(data, 'bnd')
      data = data.bnd;
    end

    if isfield(data, 'pnt')
      pos          = data.pnt;
      hastransform = 0;
      haspos       = 1;
      if isfield(data, 'tri')
        tri    = data.tri;
        hastri = 1;
      else
        hastri = 0;
      end
    else
      error('no geometrical data seems to be present in the input');
    end  
  end % switch dtype{k}
  
  if ~isfield(data, 'coordinates'), 
    data.coordinates = [];
  end
  [labelx, labely, labelz] = xyz2label(data.coordinates);
 
  % plot the geometrical object 
  
  % plot the mesh
  if haspos && hastri
    mesh.pnt = pos;
    mesh.tri = tri;
    ft_plot_mesh(mesh, 'edgecolor','none', 'facecolor', [0.6 0.8 0.6], 'facealpha', 0.6);
    camlight;
  end

  % plot the points
  if haspos && ~hastri
    shape.pnt = pos;
    ft_plot_headshape(shape);
    camlight 
  end
 
  % plot 3 slices 
  if hastransform
    ft_plot_ortho(data.anatomy, 'transform', data.transform, 'resolution', 1, 'style', 'intersect'); 
    axis vis3d
    view([110 36]);
  end

end

% get the xyz-axes
xdat  = [-axmax 0 0; axmax 0 0];  
ydat  = [0 -axmax 0; 0 axmax 0];  
zdat  = [0 0 -axmax; 0 0 axmax];  

% plot axes
hl = line(xdat, ydat, zdat);
set(hl(1), 'linewidth', 2, 'color', 'r');
set(hl(2), 'linewidth', 2, 'color', 'g');
set(hl(3), 'linewidth', 2, 'color', 'b');

% create the ball at the origin
[O.pnt, O.tri] = icosahedron42;
O.pnt = O.pnt.*rbol;
ft_plot_mesh(O, 'edgecolor', 'none');

hx = text(xdat(1,1),ydat(1,1),zdat(1,1),labelx{1},'color','y','fontsize',15,'linewidth',2);
hy = text(xdat(1,2),ydat(1,2),zdat(1,2),labely{1},'color','y','fontsize',15,'linewidth',2);
hz = text(xdat(1,3),ydat(1,3),zdat(1,3),labelz{1},'color','y','fontsize',15,'linewidth',2);
hx = text(xdat(2,1),ydat(2,1),zdat(2,1),labelx{2},'color','y','fontsize',15,'linewidth',2);
hy = text(xdat(2,2),ydat(2,2),zdat(2,2),labely{2},'color','y','fontsize',15,'linewidth',2);
hz = text(xdat(2,3),ydat(2,3),zdat(2,3),labelz{2},'color','y','fontsize',15,'linewidth',2);

% print some feedback to the user
fprintf('you can now check whether the axes and the origin are where you expect them to be. you can rotate the object, by switching on the rotation utility in the figure panel\n');

% interactively determine orientation
str{1} = input('In which direction is the positive X-axis pointing?\nr(ight),l(eft),a(nterior),p(osterior),s(uperior),i(nferior)', 's');
str{2} = input('In which direction is the positive Y-axis pointing?\nr(ight),l(eft),a(nterior),p(osterior),s(uperior),i(nferior)', 's');
str{3} = input('In which direction is the positive Z-axis pointing?\nr(ight),l(eft),a(nterior),p(osterior),s(uperior),i(nferior)', 's');

% create output
for k = 1:ndata
  varargout{k} = varargin{k};
  varargout{k}.coordinates = upper(cat(2, str{:}));
end


%-------------------------------------------------
%function to go from aplrsi to better interpretable format
function [labelx, labely, labelz] = xyz2label(str)

if ~isempty(str) && ~strcmp(str, 'unknown'),
  strx = tokenize(str, '_');
 
  switch strx{1}
  case 'ras'
    labelx = {'left'      'right'   };
    labely = {'posterior' 'anterior'};
    labelz = {'inferior'  'superior'};
  case 'als'
    labelx = {'posterior' 'anterior'};
    labely = {'right'     'left'};
    labelz = {'inferior'  'superior'};
  otherwise
  end  

else
  labelx = {'-X' 'X'};
  labely = {'-Y' 'Y'};
  labelz = {'-Z' 'Z'};
end

