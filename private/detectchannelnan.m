function [badchannels] = detectchannelnan(cfg,data)

% DETECTCHANNELNAN detects whether one or more channels contain only NaNs
% and if so adds them to cfg.badchannels
%
% Copyright (C) 2003, Robert Oostenveld
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

badchannels      = ft_getopt(cfg, 'badchannel',     {});

if any(isnan(data.trial{1}(:,1))) % check if any NaNs in first sample of first trial
    fprintf('Found channels starting with NaNs. Will be added to cfg.badchannels\n');
    badChan = ft_channelselection(find(isnan(data.trial{1}(:,1))),data.label);
    for k = 1:length(badChan)
        fprintf('Channel %s contains NaNs, added to bad channels\n',string(badChan(k)));
    end
    badchannels = cat(1, badchannels, badChan);
end
