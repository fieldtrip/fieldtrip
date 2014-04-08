function [name1, value1] = splitstruct(name0, value0)

% SPLITSTRUCT splits a structure into names and values
%
% See also PRINTSTRUCT

name1  = cell(0,1);
value1 = cell(0,1);

if isstruct(value0)
  fn = fieldnames(value0);
  
  if numel(value0)==1
    for i=1:length(fn)
      [name2, value2] = splitstruct([name0 '.' fn{i}], value0.(fn{i}));
      name1  = cat(1, name1, name2);
      value1 = cat(1, value1, value2);
    end
    clear name2 value2
  else
    siz = size(value0);
    dim = numel(siz);
    if dim>2
      error('not yet supported');
    end
    for i1=1:siz(1)
      for i2=1:siz(2)
        for i=1:length(fn)
          name2 = sprintf('%s(%d,%d)', name0, i1, i2);
          [name2, value2] = splitstruct([name2 '.' fn{i}], value0(i1,i2).(fn{i}));
          name1  = cat(1, name1, name2);
          value1 = cat(1, value1, value2);
        end
        clear name2 value2 name2
      end % i2
    end % i1
  end
  
elseif iscell(value0)
  
  siz = size(value0);
  dim = numel(siz);
  if dim>2
    error('not yet supported');
  end
  for i1=1:siz(1)
    for i2=1:siz(2)
      name2 = sprintf('%s{%d,%d}', name0, i1, i2);
      [name2, value2] = splitstruct(name2, value0{i1,i2});
      name1  = cat(1, name1, name2);
      value1 = cat(1, value1, value2);
      clear name2 value2 name2
    end % i2
  end % i1
  
else
  name1{end+1,1} = name0;
  value1{end+1,1} = value0;
end
