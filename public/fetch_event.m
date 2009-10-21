function [event]=fetch_event(data)

% FETCH_EVENT mimics the behaviour of READ_EVENT, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   [event] = fetch_event(data)
%
% See also READ_EVENT, FETCH_HEADER, FETCH_DATA

% Copyright (C) 2008, Esther Meeuwissen
%
% $Log: fetch_event.m,v $
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.2  2008/09/29 21:12:39  roboos
% cleaned up the code from Esther, added copyrights, updated documentation
%

% check whether input is data
data = checkdata(data, 'datatype', 'raw');

% locate the event structure
event = findcfg(data.cfg, 'event');

