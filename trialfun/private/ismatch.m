function s = ismatch(x, y)

% ISMATCH returns true if x is a member of array y, regardless of the class
% of x and y, if y is a string, or a cell-array of strings, it can contain
% the wildcard '*'

if isempty(x) || isempty(y)
  s = false;
  
elseif ischar(x) && ischar(y)
  y = sprintf('%s%s%s', '^', regexptranslate('wildcard',y), '$');
  s = ~isempty(regexp(x, y, 'once'));
  
elseif isnumeric(x) && isnumeric(y)
  s = ismember(x, y);
  
elseif ischar(x) && iscellstr(y)
  y = y(strcmp(class(x), cellfun(@class, y, 'UniformOutput', false)));
  s = ismember(x, y);
  % one or more of the elements in y can contain a wildcard, only proceed if s=false
  if ~s && any(contains(y, '*'))
    y = y(contains(y, '*'));
    for i = 1:numel(y)
      tmpy = sprintf('%s%s%s', '^', regexptranslate('wildcard',y{i}), '$');
      s = ~isempty(regexp(x, tmpy, 'once'));
      if s, return; end
    end
  end
  
elseif isnumeric(x) && iscell(y)
  % this works if y contains both numbers and strings
  s = false;
  for i=1:numel(y)
    s = s || ismember(x, y{i});
  end
  
elseif ischar(x) && iscell(y)
  % this works if y contains both numbers and strings
  s = false;
  for i=1:numel(y)
    s = s || strcmp(x, y{i});
  end
  
else
  s = false;
end
