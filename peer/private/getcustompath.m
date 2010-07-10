%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that determines the path, excluding all Matlab toolboxes
% the directories and the path on a windows computer look different than on unix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, d] = getcustompath

% these are for faster processing on subsequent calls
persistent previous_path previous_argout

if isequal(path, previous_path)
  % don't do the processing again, but return the previous values from cache
  p = previous_argout{1};
  d = previous_argout{2};
  return
end

if ispc
  p = tokenize(path, ';');
else
  p = tokenize(path, ':');
end
% remove the matlab specific directories
if ispc
  s = false(size(p));
  for i=1:length(p)
    s(i) = ~strncmp(p{i}, matlabroot, length(matlabroot));
  end
else
  s = cellfun(@isempty, regexp(p, ['^' matlabroot]));
end
d = p(~s);
p = p( s);
% remove the directory containing the peer code, the slave should use its own
p = setdiff(p, fileparts(mfilename('fullpath')));
% concatenate the path, using the platform specific seperator
if ispc
  p = sprintf('%s;', p{:});
  d = sprintf('%s;', d{:});
else
  p = sprintf('%s:', p{:});
  d = sprintf('%s:', d{:});
end
p = p(1:end-1); % remove the last separator
d = d(1:end-1); % remove the last separator

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_path   = path;
previous_argout = {p, d};

