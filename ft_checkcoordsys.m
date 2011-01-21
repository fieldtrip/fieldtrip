function [data] = ft_checkcoordsys(cfg, data)

% FT_CHECKCOORDSYS allows the user to perform  a check on the coordinatesystem
% and units in which the input object is expressed.
%
% Use as
%   [output] = ft_checkcoordsys(cfg, data)
%
% See also FT_VOLUMEREALIGN, FT_VOLUMERESLICE

if ~isfield(cfg, 'interactive'), cfg.interactive = 'yes'; end

dointeractive = strcmp(cfg.interactive, 'yes');
doplot        = 1;

dtype = datatype(data);
data  = ft_convert_units(data);
unit  = data.unit;

% determine the size of the "unit" sphere in the origin and the length of the axes
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
    error('unknown units (%s)', unit);
end

% the plotting style depends on the data content
haspos        = false;
hastri        = false;
hastransform  = false;

% do some checks on the geometrical object
switch dtype
  case 'volume'
    % contains transform
    hastransform = true;
  case 'source'
    % contains pos
    haspos       = true;
    hastri       = false; % FIXME may need to be supported in the future
  case 'dip'
    % contains pos
    haspos       = true;
    hastri       = false; % FIXME may need to be supported in the future
  case {'freq', 'raw', 'timelock', 'spike', 'mvar', 'freqmvar'}
    % these may contain geometrical objects in grad or elec fields
    if isfield(data, 'grad')
      data = data.grad;
      haspos = true;
    elseif isfield(data, 'elec')
      data = data.elec;
      haspos = true;
    end
  case 'unknown'
    % other geometrical objects typically contain geometric information in the bnd.pnt or pnt field
    if isfield(data, 'bnd')
      data = data.bnd;
    end
    
    if isfield(data, 'pnt')
      haspos       = true;
    end
    if isfield(data, 'tri')
      hastri = true;
    end
end % switch dtype{k}

if haspos
  if isfield(data, 'pnt')
    pos = data.pnt;
  else isfield(data, 'pos')
    pos = data.pos;
  end
end

if hastri
  tri = data.tri;
end

if ~isfield(data, 'coordsys'),
  data.coordsys = [];
end

if doplot,
  % create the labels that are to be plotted along the axes
  [labelx, labely, labelz] = xyz2label(data.coordsys);
  
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
  
  % get the xyz-axes
  xdat  = [-axmax 0 0; axmax 0 0];
  ydat  = [0 -axmax 0; 0 axmax 0];
  zdat  = [0 0 -axmax; 0 0 axmax];
  
  % get the xyz-axes dotted
  xdatdot = (-axmax:(axmax/15):axmax);
  xdatdot = xdatdot(1:floor(numel(xdatdot)/2)*2);
  xdatdot = reshape(xdatdot, [2 numel(xdatdot)/2]);
  n       = size(xdatdot,2);
  ydatdot = [zeros(2,n) xdatdot zeros(2,n)];
  zdatdot = [zeros(2,2*n) xdatdot];
  xdatdot = [xdatdot zeros(2,2*n)];
  
  % plot axes
  hl = line(xdat, ydat, zdat);
  set(hl(1), 'linewidth', 1, 'color', 'r');
  set(hl(2), 'linewidth', 1, 'color', 'g');
  set(hl(3), 'linewidth', 1, 'color', 'b');
  hld = line(xdatdot, ydatdot, zdatdot);
  for k = 1:n
    set(hld(k    ), 'linewidth', 3, 'color', 'r');
    set(hld(k+n*1), 'linewidth', 3, 'color', 'g');
    set(hld(k+n*2), 'linewidth', 3, 'color', 'b');
  end
  
  % create the ball at the origin
  [O.pnt, O.tri] = icosahedron42;
  O.pnt = O.pnt.*rbol;
  ft_plot_mesh(O, 'edgecolor', 'none');
  
  % add the labels to the axis
  text(xdat(1,1),ydat(1,1),zdat(1,1),labelx{1},'color','y','fontsize',15,'linewidth',2);
  text(xdat(1,2),ydat(1,2),zdat(1,2),labely{1},'color','y','fontsize',15,'linewidth',2);
  text(xdat(1,3),ydat(1,3),zdat(1,3),labelz{1},'color','y','fontsize',15,'linewidth',2);
  text(xdat(2,1),ydat(2,1),zdat(2,1),labelx{2},'color','y','fontsize',15,'linewidth',2);
  text(xdat(2,2),ydat(2,2),zdat(2,2),labely{2},'color','y','fontsize',15,'linewidth',2);
  text(xdat(2,3),ydat(2,3),zdat(2,3),labelz{2},'color','y','fontsize',15,'linewidth',2);
end

if dointeractive,
  % print some feedback to the user
  fprintf('You should now check whether the axes and the origin are consistent with what you expect them to be\n');
  fprintf('Switch on the rotation utility in the figure panel to see the figure from all angles\n');
  fprintf('For each of the following questions press one of the buttons: r(ight),l(eft),a(nterior),p(osterior),s(uperior),i(nferior)\n');
  
  % interactively determine orientation
  str = {'', '', ''};
  while ~ismember(str{1}, {'r', 'l', 'a', 'p', 's', 'i'})
    str{1} = input('In which direction is the positive X-axis pointing? ', 's');
  end
  while ~ismember(str{2}, {'r', 'l', 'a', 'p', 's', 'i'})
    str{2} = input('In which direction is the positive Y-axis pointing? ', 's');
  end
  while ~ismember(str{3}, {'r', 'l', 'a', 'p', 's', 'i'})
    str{3} = input('In which direction is the positive Z-axis pointing? ', 's');
  end
  
  % interactively determine origin
  str2 = input('Does the origin of the coordinate system look as any of the following?\na(nterior commissure),i(nterauricular),n(ot a landmark)', 's');
  
  % create output
  str = upper(cat(2, str{:}));
  switch str2
    case 'a'
      if strcmp(str, 'RAS')
        str = [str,'_TAL'];
      else
      end
    case 'i'
      if strcmp(str, 'ALS')
        str = [str,'_CTF'];
      elseif strcmp(str, 'RAS')
        str = [str,'_Neuromag'];
      end
    otherwise
  end
  
  data.coordsys = str;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to go from aplrsi to better interpretable format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [labelx, labely, labelz] = xyz2label(str)

if ~isempty(str) && ~strcmp(str, 'unknown'),
  strx = tokenize(str, '_');
  
  switch lower(strx{1})
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
  labelx = {'-X' '+X'};
  labely = {'-Y' '+Y'};
  labelz = {'-Z' '+Z'};
end
