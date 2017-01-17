function table = struct2table(s)

% STRUCT2TABLE converts a struct-array to a string that contains a pretty-print table
%
% Example
%   s(1).a = 1
%   s(1).b = 2
%   s(2).a = 3
%   s(2).b = 4
%   disp(struct2table(s))

fn = fieldnames(s);
fl = zeros(size(fn));

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
  fl(i) = max(cellfun(@length, {s.(fn{i})}));
  fl(i) = max(fl(i), length(fn{i}));
end

header = '|';
for i=1:numel(fn)
  header = [header pad(fn{i},fl(i)+1,'left',' '), ' |'];
end


divider = repmat('-', size(header));

% divider = header;
% divider(divider~='|') = '-';

line = cell(numel(s),1);
for i=1:numel(s)
  line{i} = '|';
  for j=1:numel(fn)
    line{i} = [line{i} pad(s(i).(fn{j}),fl(j)+1,'left',' '), ' |'];
  end
end

table = cat(1, {header, divider}', line);

