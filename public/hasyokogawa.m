function [version] = hasyokogawa(desired)

% HASYOKOGAWA tests whether the data input toolbox for MEG systems by
% Yokogawa (www.yokogawa.com, designed by KIT/EagleTechnology) is
% installed. Only the newest version of the toolbox is accepted.
%
% Use as
%   [string]  = hasyokogawa;
% which returns a string describing the toolbox version, e.g. "12bitBeta3",
% "16bitBeta3", or "16bitBeta6". An empty string is returned if the toolbox
% is not installed. The string "unknown" is returned if it is installed but
% the version is unknown.
%
% Alternatively you can use it as
%   [boolean] = hasyokogawa(desired);
% where desired is a string with the desired version.
%
% See also READ_YOKOGAWA_HEADER, READ_YOKOGAWA_DATA, READ_YOKOGAWA_EVENT,
% YOKOGAWA2GRAD

% Copyright (C) 2010, Tilmann Sander-Thoemmes
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

% return empty if not present
version = [];

try
  warning('off', 'MATLAB:pfileOlderThanMfile');
  % Call some functions with input argument "Inf": If
  % the functions are present they return their revision number.
  % Call first GetMeg160ADbitInfoM as this is not present in
  % the 12bit library, in case of error the "catch" part will take over.
  % The code below is intentionally very literal for easy of reading.
  res = textscan(evalc('GetMeg160ADbitInfoM(Inf);'),'%s %s %c %s %s %d');
  rev_ADbitInfoM = res{6};
  res = textscan(evalc('GetMeg160ChannelInfoM(Inf);'),'%s %s %c %s %s %d');
  rev_ChannelInfoM = res{6};
  res = textscan(evalc('GetMeg160AmpGainM(Inf);'),'%s %s %c %s %s %d');
  rev_AmpGainM = res{6};
  res = textscan(evalc('GetMeg160MatchingInfoM(Inf);'),'%s %s %c %s %s %d');
  rev_MatchingInfoM = res{6};
  if [0 2 1 5] == [rev_ADbitInfoM rev_ChannelInfoM rev_AmpGainM rev_MatchingInfoM]
    version='16bitBeta3';
  elseif [0 2 2 5] == [rev_ADbitInfoM rev_ChannelInfoM rev_AmpGainM rev_MatchingInfoM]
    version='16bitBeta6';
  else
    warning('Yokogawa toolbox is installed, but the version cannot be determined.');
    version = 'unknown';
  end
  if nargin>0
    version = strcmpi(version, desired);
    if ~version
        warning('The required version of the Yokogawa input toolbox (%s) is not installed.', desired);
    end
  end
  warning('on', 'MATLAB:pfileOlderThanMfile');
catch
  warning('on', 'MATLAB:pfileOlderThanMfile');
  m = lasterror;
  m.identifier;
  if strcmp(m.identifier, 'MATLAB:UndefinedFunction') || strcmp(m.identifier, 'MATLAB:FileIO:InvalidFid')
    if (exist('GetMeg160ChannelInfoM') && exist('GetMeg160ContinuousRawDataM'));
      version = '12bitBeta3';
    else
      version = 'unknown';
   end
  end
  if nargin>0
    version = 0; % logical output
  end
  
end
