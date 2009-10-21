function [selected] = select_point(pos, varargin)

% SELECT_POINT helper function for selecting a one or multiple points
% in the current figure using the mouse.
%
% Use as
%   [selected] = select_point(pos, ...)
%
% It returns a list of the [x y] coordinates of the selected points.
%
% Optional input arguments should come in key-value pairs and can include
%   'multiple'    true/false, make multiple selections, pressing "q" on the keyboard finalizes the selection (default = false)
%   'nearest'     true/false (default = true)
%
% Example use
%   pos = randn(10,2);
%   figure
%   plot(pos(:,1), pos(:,2), '.')
%   select_point(pos)

% $Log: select_point.m,v $
% Revision 1.2  2009/06/30 11:46:15  roboos
% fixed docu
%
% Revision 1.1  2009/06/30 11:44:34  roboos
% renamed select_pointd into select_point for consistency with plot_topo
%
% Revision 1.7  2009/06/22 12:33:11  crimic
% minor change
%
% Revision 1.6  2009/06/16 08:17:40  crimic
% added check on input
%
% Revision 1.5  2009/06/15 15:46:45  roboos
% first implementation of point3d, multiple changes to point2d, still some work to be done to make them consistent
%
% Revision 1.4  2009/06/15 13:43:27  roboos
% reimplemented from scratch
%
% Revision 1.3  2009/06/04 10:51:09  roboos
% only whitespace
%
% Revision 1.2  2009/06/03 11:23:46  crimic
% first implementation
%


% get optional input arguments
nearest  = keyval('nearest', varargin); if isempty(nearest), nearest = true; end
multiple = keyval('multiple', varargin); if isempty(multiple), multiple = false; end

% ensure that it is boolean
nearest = istrue(nearest);
multiple = istrue(multiple);

if multiple
  fprintf('select multiple points by clicking in the figure, press "q" if you are done\n');
end

x = [];
y = [];
done = false;

selected = zeros(0,3);

% ensure that "q" is not the current character, which happens if you reuse the same figure
set(gcf, 'CurrentCharacter', 'x')

while ~done
  k     = waitforbuttonpress;
  point = get(gca,'CurrentPoint');     % button down detected
  key   = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
  if strcmp(key, 'q')
    % we are done with the clicking
    done = true;
  else
    % add the current point
    x(end+1) = point(1,1);
    y(end+1) = point(1,2);
  end

  if ~multiple
    done = true;
  end
end

if nearest && ~isempty(pos)
  % determine the points that are the nearest to the displayed points
  selected = [];

  % compute the distance between the points to get an estimate of the tolerance
  dp = dist(pos');
  dp = triu(dp, 1);
  dp = dp(:);
  dp = dp(dp>0);
  % allow for some tolerance in the clicking
  dp = median(dp);
  tolerance = 0.3*dp;

  for i=1:length(x)
    % compute the distance between the clicked position and all points
    dx = pos(:,1) - x(i);
    dy = pos(:,2) - y(i);
    dd = sqrt(dx.^2 + dy.^2);
    [d, i] = min(dd);
    if d<tolerance
      selected(end+1,:) = pos(i,:);
    end
  end

else
  selected = [x(:) y(:)];
end % if nearest
