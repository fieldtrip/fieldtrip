function offset = time2offset(time, fsample);

% TIME2OFFSET converts a time-axis of a trial into the offset in samples
% according to the definition from DEFINETRIAL
%
% Use as
%   [offset] = time2offset(time, fsample)
%
% The trialdefinition "trl" is an Nx3 matrix. The first column contains
% the sample-indices of the begin of the trial relative to the begin
% of the raw data , the second column contains the sample_indices of
% the end of the trials, and the third column contains the offset of
% the trigger with respect to the trial. An offset of 0 means that
% the first sample of the trial corresponds to the trigger. A positive
% offset indicates that the first sample is later than the triger, a
% negative offset indicates a trial beginning before the trigger.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: time2offset.m,v $
% Revision 1.1  2008/11/20 13:48:47  roboos
% moved from private to public
%
% Revision 1.2  2005/08/05 09:18:21  roboos
% round the offset to the nearest integer
% added copyright and cvs log
%

offset = round(time(1)*fsample);
