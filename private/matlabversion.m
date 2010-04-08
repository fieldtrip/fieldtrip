function [v] = matlabversion;

% MATLABVERSION returns the Matlab version as a number
%
% Use as
%  v = matlabversion;
%
% An example for using this function is
%
%  if matlabversion<7
%    % do something specific to support an old matlab version
%  else
%    % do the default operation
%  end
%
% See also VERSION, VER

% Copyright (C) 2006, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

s = ver('matlab');
v = s.Version;

if ischar(v) 
  % try converting to a number
  n = str2num(v);
  if isempty(n)
    switch v
    case '6.5.1'
      n = 6.5; % this is accurate enough
    case '7.0.1'
      n = 7.0; % this is accurate enough
    case '7.0.4'
      n = 7.0; % this is accurate enough
    otherwise  
      warning('cannot convert matlab version into a number');
      v = v;
    end
  end
  v = n;
end

