function mutexunlock(lockfile)

% MUTEXUNLOCK removes a lockfile
%
% Use as
%   mutexunlock(lockfile)
%
% See also MUTEXLOCK and http://en.wikipedia.org/wiki/Mutual_exclusion

if isempty(lockfile)
  % no lockfile was specified, return immediately
  return
end

delete(lockfile);

