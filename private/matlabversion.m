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
% $Log: matlabversion.m,v $
% Revision 1.3  2007/07/07 12:10:02  roboos
% detect matlab version 7.0.1 and 7.0.4 and treat as 7.0
%
% Revision 1.2  2006/06/08 07:08:38  roboos
% convert string into number
%
% Revision 1.1  2006/05/29 07:58:25  roboos
% new helper function
%

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

