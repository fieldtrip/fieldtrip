function str = printand(strs, quote)

% PRINTAND prints a single or mutiple strings as "x1, x2, x3 and x4". If there is
% only one string, that string is returned without additional formatting.
%
% See also PRINTOR

if nargin<2
  quote = false;
end

if quote
  if numel(strs)>2
    str = sprintf('''%s'', ', strs{1:(end-2)});
    str = sprintf('%s''%s'' and ''%s''', str, strs{end-1}, strs{end});
  elseif numel(strs)==2
    str = sprintf('''%s'' and ''%s''', strs{end-1}, strs{end});
  else
    str = strs{1};
  end
else
  if numel(strs)>1
    str = sprintf('%s, ', strs{1:(end-2)});
    str = sprintf('%s%s and %s', str, strs{end-1}, strs{end});
  else
    str = strs{1};
  end
end