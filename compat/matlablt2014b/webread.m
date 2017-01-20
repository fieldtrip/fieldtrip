function s = webread(varargin)

options = {};
url = varargin{1};
s = urlread(url, 'get', options);

if isequal(s(1:2), '[{')
  % assume that it is json formatted
  split = regexp(s, '},{');
  b = [2 split+1];
  e = [split length(s)-1];
  c = cell(size(b));
  for i=1:numel(c)
    c{i} = json2struct(s(b(i):e(i)));
  end
  if numel(c)>1
    % return a cell-array with structs
    s = c;
  elseif numel(c)==1
    % return a single struct
    s = c{1}; 
  else
    % return the string as it is
  end
end
