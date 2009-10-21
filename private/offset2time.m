function time = offset2time(offset, fsample, nsamples);

% OFFSET2TIME converts the offset of a trial definition into a time-axis
% according to the definition from DEFINETRIAL
%
% Use as
%   [time] = offset2time(offset, fsample, nsamples)
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
% $Log: offset2time.m,v $
% Revision 1.2  2009/02/06 10:10:47  roboos
% convert nsamples and offset to double precision
%
% Revision 1.1  2008/11/20 13:48:47  roboos
% moved from private to public
%
% Revision 1.2  2005/08/05 09:19:00  roboos
% added copyright and cvs log
%

% ensure that these are not integers
offset   = double(offset);
nsamples = double(nsamples);

time = (offset + (0:(nsamples-1)))/fsample;
