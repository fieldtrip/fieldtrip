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

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information


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
