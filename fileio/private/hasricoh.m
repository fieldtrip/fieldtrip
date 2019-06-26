function [version] = hasricoh(desired)

%% HASRICOH tests whether the official toolbox for RICOH MEG systems by
%% Ricoh Company, Ltd. is installed or not.
%% Use as
%%   string  = hasricoh;
%% which returns a string describing the toolbox version '1.0'.
%% An empty string is returned if the toolbox is not installed.
%% The string "unknown" is returned if it is installed but
%% the version is unknown.
%%
%% Alternatively you can use it as
%%   [boolean] = hasricoh(desired);
%% where desired is a string with the desired version.
%%
%% See also READ_RICOH_HEADER, READ_RICOH_DATA, READ_RICOH_EVENT, RICOH2GRAD

if exist('getRVersion', 'file')
  res = getRVersion();
  version = res.version;
else
  % return empty if none of them is present
  version = [];
end

if nargin>0
  % return a true/false value
  if isempty(version)
    version = false;
  else
    version = strcmpi(version, desired);
  end
end
