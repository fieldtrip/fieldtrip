function str = printor(strs)

% PRINTOR prints a single or mutiple strings as "x1, x2, x3 or x4". If there is
% only one string, that string is returned without additional formatting.
%
% See also PRINTOR

if numel(strs)>1
  str = sprintf('%s, ', strs{1:(end-2)});
  str = sprintf('%s%s or %s', str, strs{end-1}, strs{end});
else
  str = strs{1};
end
