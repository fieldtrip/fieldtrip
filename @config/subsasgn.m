function y = subsasgn(x, index, val)

% SUBSASGN Assign a new value to a specified field in a config objects and increment its assignment counter.

if length(index)==1
  switch index.type
    case '.'
      y = set(x, index.subs, val);
    case '{}'
      error('Cell contents reference from a non-cell array object.');
    case '()'
      error('Index exceeds matrix dimensions.');
    otherwise
      error('Incorrect contents reference');
  end
elseif length(index)>1 && strcmp(index(1).type, '()')
  if ~isfield(x, index(2).subs)
    % an empty field should be added to all elements of the array
    for i=1:numel(x)
      y(i) = subsasgn(x(i), index(2:end), []);
    end
  else
    % the field is already present in the array
    y = x;
  end
  % the value of the field should only be changed for the specific element of the array
  y(index(1).subs{1}) = subsasgn(y(index(1).subs{1}), index(2:end), val);
else
  % use recursion to find the subfield that is being indexed
  if ~isfield(x, index(1).subs)
    % the subfield does not exists, create a new config subfield
    x  = subsasgn(x, index(1), config());
  end
  y1 = subsref(x, index(1), false);
  y2 = subsasgn(y1, index(2:end), val);
  y  = set(x, index(1).subs, y2);
end

