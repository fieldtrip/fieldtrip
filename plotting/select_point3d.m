function [selected] = select_point3d(bnd, varargin)

% SELECT_POINT3D helper function for selecting one or multiple points
% on a 3D mesh using the mouse.
%
% Use as
%   [selected] = select_point3d(bnd, ...)
%
% It returns a list of the [x y z] coordinates of the selected points.
%
% Optional input arguments should come in key-value pairs and can include
%   'multiple'    true/false, make multiple selections, pressing "q" on the keyboard finalizes the selection (default = false)
%   'nearest'     true/false (default = true)
%
% Example use
%   [pnt, tri] = icosahedron162;
%   bnd.pnt = pnt;
%   bnd.tri = tri;
%   plot_mesh(bnd)
%   camlight
%   ... do something here

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get optional input arguments
nearest  = keyval('nearest', varargin); if isempty(nearest), nearest = true; end
multiple = keyval('multiple', varargin); if isempty(multiple), multiple = false; end

% ensure that it is boolean
nearest = istrue(nearest);
multiple = istrue(multiple);

% get the object handles
h = get(gca, 'children');

% select the correct objects
iscorrect = false(size(h));
for i=1:length(h)
  try
    pnt = get(h(i),'vertices');
    tri = get(h(i),'faces');
    if ~isempty(bnd) && isequal(bnd.pnt, pnt) && isequal(bnd.tri, tri)
      % it is the same object that the user has plotted before
      iscorrect(i) = true;
    elseif isempty(bnd)
      % assume that it is the same object that the user has plotted before
      iscorrect(i) = true;
    end
  end
end
h = h(iscorrect);

if isempty(h) && ~isempty(bnd)
  figure
  plot_mesh(bnd);
  camlight
  selected = select_point3d(bnd, varargin{:});
  return
end

if length(h)>1
  warning('using the first patch object in the figure');
  h = h(1);
end

selected = zeros(0,3);

done = false;
while ~done
  k = waitforbuttonpress;
  [p v vi facev facei] = select3d(h);
  key = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?

  if strcmp(key, 'q')
    % finished selecting points
    done = true;
  else
    % a new point was selected
    if nearest
      selected(end+1,:) = v;
    else
      selected(end+1,:) = p;
    end % if nearest
    fprintf('selected point at [%f %f %f]\n', selected(end,1), selected(end,2), selected(end,3));
  end

  if ~multiple
    done = true;
  end
end

