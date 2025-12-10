function [trl] = continuous_ns(cfg)

% CONTINUOUS_NS created a trial definition from a Neuroscan *.cnt file
% which subsequently can be used in the EEG/MEG framework
%
% Use as
%   [trl] = continuous_ns(cfg)
%
% where the configuration should contain 
%   cfg.trialdef.trigger  = number or list with triggers
%   cfg.trialdef.prestim  = pre-stimulus in seconds
%   cfg.trialdef.poststim = post-stimulus in seconds    
%
% See also SINGLETRIAL_NS

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

% read the header and the event table from the continuous file
tmp = read_ns_cnt(cfg.datafile, 'ldheaderonly', 1);
fprintf('found %d events in continuous neuroscan file\n', length(tmp.event.frame));
eventindx = find(ismember(tmp.event.stimtype, cfg.trialdef.trigger));
trl(:,1) = tmp.event.frame(eventindx) - round(cfg.trialdef.prestim*tmp.rate) + 1;  % begin sample
trl(:,2) = tmp.event.frame(eventindx) + round(cfg.trialdef.poststim*tmp.rate) + 1; % end sample   
trl(:,3) = -round(cfg.trialdef.prestim*tmp.rate);                                  % offset
fprintf('selected %d events based on triggercode\n', size(trl,1));

return
