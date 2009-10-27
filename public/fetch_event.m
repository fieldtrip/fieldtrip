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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% check whether input is data
data = checkdata(data, 'datatype', 'raw');

% locate the event structure
event = findcfg(data.cfg, 'event');

