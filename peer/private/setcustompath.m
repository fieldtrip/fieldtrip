function setcustompath(option_path)

% SUBFUNCTION that changes the path 

% these are for remembering the path on subsequent calls with the same input arguments
persistent previous_argin previous_path

if isequal(previous_argin, option_path) && isequal(previous_path, path)
  % no reason to change the path again
  return
end

if ~isempty(option_path)
  [p, d] = getcustompath;
  % ensure that the default path is maintained
  % and that the peer functions are maintained on the path
  path(d);
  path(option_path, path);
  path(fileparts(which('peerfeval')), path);
end

% remember the current settings for the next call
previous_argin = option_path;
previous_path  = path;

