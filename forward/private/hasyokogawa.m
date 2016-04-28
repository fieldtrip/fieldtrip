function [version] = hasyokogawa(desired)

% HASYOKOGAWA tests whether the data input toolbox for MEG systems by
% Yokogawa (www.yokogawa.com, designed by KIT/EagleTechnology) is
% installed. Only the newest version of the toolbox is accepted.
%
% Use as
%   string  = hasyokogawa;
% which returns a string describing the toolbox version, e.g. "12bitBeta3",
% "16bitBeta3", or "16bitBeta6" for preliminary versions, or '1.4' for the
% official Yokogawa MEG Reader Toolbox. An empty string is returned if the toolbox
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
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

ws = warning('off', 'MATLAB:pfileOlderThanMfile');

% there are a few versions of the old preliminary implementation, such as
% 12bitBeta3, 16bitBeta3 and 16bitBeta6. In 2011 a completely new
% implementation was officially released, which contains functions with
% other names. At the time of writing this, the new implementation is
% version 1.4.

if exist('getYkgwVersion', 'file')
  res = getYkgwVersion();
  version = res.version;

elseif exist('GetMeg160ADbitInfoM', 'file') || exist('GetMeg160ChannelInfoM', 'file') || exist('GetMeg160ContinuousRawDataM', 'file')
  % start with unknown, try to refine the version
  version = 'unknown';

  try
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
      warning('The version of the installed Yokogawa toolbox cannot be determined.');
    end
  catch
    m = lasterror;
    m.identifier;
    if strcmp(m.identifier, 'MATLAB:UndefinedFunction') || strcmp(m.identifier, 'MATLAB:FileIO:InvalidFid')
      if (exist('GetMeg160ChannelInfoM', 'file') && exist('GetMeg160ContinuousRawDataM', 'file'));
        version = '12bitBeta3';
      end
    end
  end

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

% revert to the original warning state
warning(ws);
