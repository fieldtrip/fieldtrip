function tf = endsWith(s, pattern, str, boolean)

%ENDSWITH True if text starts with pattern.
%   TF = endsWith(S,PATTERN) returns true if any element of string array S
%   starts with PATTERN. TF is the same size as S.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2016b.

if ~ischar(s) && ~iscellstr(s)
  error('the input should be either a char-array or a cell array with chars');
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

% the final comparison is done on the start of the string, hence flip all of them
pattern = fliplr(pattern);
if iscell(s)
  s = cellfun(@fliplr, s);
else
  s = fliplr(s);
end

% compare the start of the string
if boolean
  tf = strncmpi(s, pattern, numel(pattern));
else
  tf = strncmp(s, pattern, numel(pattern));
end


