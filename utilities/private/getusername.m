function user = getusername()

% USERNAME
%
% Use as
%   str = username;

% this is to speed up subsequent calls
persistent previous_argout
if ~isempty(previous_argout)
  user = previous_argout;
  return
end

if (ispc())
  user = getenv('UserName');
else
  user = getenv('USER');
end

if (isempty(user))
  user = 'unknown';
end

% remember for subsequent calls
previous_argout = user;

