function s = json2struct(j)

j = j(2:end-1); % remove the { and }
t = tokenize(j, ',');

s = struct();
for i=1:numel(t)
  kv = tokenize(t{i}, ':');
  key = kv{1}(2:end-1); % remove the " and "
  val = kv{2};
  if any(val(1)=='0123456789') && ~isnan(str2double(val))
    val = str2double(val);
  elseif isequal(val, 'true')
    val = true;
  elseif isequal(val, 'false')
    val = false;
  else
    val = val(2:end-1); % remove the " and "
  end
  key = fixname(key);
  s.(key) = val;
end


