function table = struct2table(s)

% STRUCT2TABLE converts a struct-array to a cell-array of strings that represents a table
%
% Example
%   s(1).a = 1
%   s(1).b = 2
%   s(2).a = 3
%   s(2).b = 4
%   disp(struct2table(s))

fn = fieldnames(s);
colwidth = zeros(size(fn));

% first ensure that all fields are strings
for i=1:numel(s)
  for j=1:numel(fn)
    val = s(i).(fn{j});
    switch class(val)
      case 'char'
        % nothing to be done
      case {'single', 'double'}
        val = num2str(val);
      case 'logical'
        if val
          val = 'true';
        else
          val = 'false';
        end
      otherwise
        error('not yet implemented for class "%s"', class(val));
    end
    s(i).(fn{j}) = val;
  end
end

for i=1:numel(fn)
  colwidth(i) = max(cellfun(@length, {s.(fn{i})}));
end

header = '|';
for i=1:numel(fn)
  colname = fn{i};
  if numel(colname)>1 && colname(1)=='X' && colname(end)=='X'
    % the name of the column is base64 encoded
    colname = fixname(fn{i});
  end
  % update the width of the column to the header
  colwidth(i) = max(colwidth(i), length(colname));
  header = [header pad(colname, colwidth(i)+1, 'left', ' '), ' |'];
end

divider = repmat('-', size(header));

% divider = header;
% divider(divider~='|') = '-';

line = cell(numel(s),1);
for i=1:numel(s)
  line{i} = '|';
  for j=1:numel(fn)
    line{i} = [line{i} pad(s(i).(fn{j}),colwidth(j)+1,'left',' '), ' |'];
  end
end

table = cat(1, header, divider, line);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
% it is possible that the column name is not valid as field name
% in that case it is base64-encoded and padded with 'X'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = fixname(str)
str = str(2:end-1);   % it starts and ends with 'X'
str(str=='_') = '=';  % the '=' sign has been replaced with '_'
str = char(base64decode(str));
