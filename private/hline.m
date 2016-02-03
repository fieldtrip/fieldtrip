function hline(y, varargin)

% HLINE plot a horizontal line in the current graph

abc = axis;
y = [y y];
x = abc([1 2]);
if length(varargin)==1
  varargin = {'color', varargin{1}};
end
h = line(x, y);
set(h, varargin{:});
