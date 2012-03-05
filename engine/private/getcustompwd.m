%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that determines the present working directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = getcustompwd

% these are for faster processing on subsequent calls
persistent previous_pwd previous_argout

if isequal(pwd, previous_pwd)
  % don't do the processing again, but return the previous values from cache
  p = previous_argout;
  return
end

% don't use the present directory if it contains the peer code
% it will confuse the slave with a potentially different mex file
if strcmp(pwd, fileparts(mfilename('fullpath')))
  warning_once('the peer slave will not change directory to %s', pwd);    
  p = [];
else
  p = pwd;
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
previous_pwd    = pwd;
previous_argout = p;

