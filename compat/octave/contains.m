function tf = contains(s, pattern, str, boolean)

%CONTAINS True if text contains a pattern.
%   TF = contains(S,PATTERN) returns true if any element of string array S
%   contains PATTERN. TF is the same size as S.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2016b.

if ~ischar(s) && ~iscellstr(s)
  error('the input should be either a char-array or a cell-array with chars');
end

if nargin<4
  boolean = false;
end
if nargin<3
  str = 'IgnoreCase';
end
if ~strcmpi(str, 'ignorecase')
  error('incorrect third input argument, can only be ''IgnoreCase''');
end
if ~islogical(boolean)
  error('fourth input argument should be a logical scalar');
end

if ~iscellstr(s)
  s = {s};
end

if boolean
 s       = lower(s);
 pattern = lower(pattern);
end

tf = ~cellfun(@isempty, strfind(s, pattern));

