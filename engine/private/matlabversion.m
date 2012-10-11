function [inInterval] = matlabversion(min, max)

% MATLABVERSION checks if the current matlab version is within the interval
% specified by min and max.
%
% Use, e.g., as:
%  if matlabversion(7.0, 7.9)
%    % do something
%  end
%
% Both strings and numbers, as well as infinities, are supported, eg.:
%  matlabversion(7.1, 7.9)    % is version between 7.1 and 7.9?
%  matlabversion(6, '7.10')   % is version between 6 and 7.10? (note: '7.10', not 7.10)
%  matlabversion(-Inf, 7.6)   % is version <= 7.6?
%  matlabversion('2008b', '2010a')
%  matlabversion('2008b', Inf)
%  matlabversion('2009b')     % exactly 2009b
%  etc.
%
% See also VERSION, VER, VERLESSTHAN

% Copyright (C) 2006, Robert Oostenveld
% Copyright (C) 2010, Eelke Spaak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% this does not change over subsequent calls, making it persistent speeds it up
persistent curVer

if nargin<2
  max = min;
end

if isempty(curVer)
  curVer = version();
end

if ((ischar(min) && isempty(str2num(min))) || (ischar(max) && isempty(str2num(max))))
  % perform comparison with respect to release string
  
  ind = strfind(curVer, '(R');
  [year, ab] = parseMatlabRelease(curVer((ind + 2):(numel(curVer) - 1)));
  
  [minY, minAb] = parseMatlabRelease(min);
  [maxY, maxAb] = parseMatlabRelease(max);
  
  inInterval = orderedComparison(minY, minAb, maxY, maxAb, year, ab);
else % perform comparison with respect to version number
  [major, minor] = parseMatlabVersion(curVer);
  [minMajor, minMinor] = parseMatlabVersion(min);
  [maxMajor, maxMinor] = parseMatlabVersion(max);
  
  inInterval = orderedComparison(minMajor, minMinor, maxMajor, maxMinor, major, minor);
end

  function [year, ab] = parseMatlabRelease(str)
    if (str == Inf)
      year = Inf; ab = Inf;
    elseif (str == -Inf)
      year = -Inf; ab = -Inf;
    else
      year = str2num(str(1:4));
      ab = str(5);
    end
  end

  function [major, minor] = parseMatlabVersion(ver)
    if (ver == Inf)
      major = Inf; minor = Inf;
    elseif (ver == -Inf)
      major = -Inf; minor = -Inf;
    elseif (isnumeric(ver))
      major = floor(ver);
      minor = int8((ver - floor(ver)) * 10);
    else % ver is string (e.g. '7.10'), parse accordingly
      [major, rest] = strtok(ver, '.');
      major = str2num(major);
      minor = str2num(strtok(rest, '.'));
    end
  end

% checks if testA is in interval (lowerA,upperA); if at edges, checks if testB is in interval (lowerB,upperB).
  function inInterval = orderedComparison(lowerA, lowerB, upperA, upperB, testA, testB)
    if (testA < lowerA || testA > upperA)
      inInterval = false;
    else
      inInterval = true;
      if (testA == lowerA)
        inInterval = inInterval && (testB >= lowerB);
      end
      
      if (testA == upperA)
        inInterval = inInterval && (testB <= upperB);
      end
    end
  end

end % function
