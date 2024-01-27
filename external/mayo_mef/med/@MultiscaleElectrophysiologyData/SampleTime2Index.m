function [sample_index, sample_yn] = SampleTime2Index(this, sample_time, options)
% MULTISCALEELECTROPHYSIOLOGYDATA.SAMPLETIME2INDEX Convert sample time to sample index
%
% Syntax:
%   [sample_index, sample_yn] = SampleTime2Index(this, sample_time)
%   [sample_index, sample_yn] = SampleTime2Index(__, st_unit)
%
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyData object
%   sample_time     - [num array] array of sample time (default unit uUTC)
%   st_unit         - [str] (optional) sample time unit: 'uUTC' (default)
%                     or 'u', 'mSec', 'Second' or 's', 'Minute' or 'm', 'Hour' or
%                     'h' and 'Day' or 'd'.
%
% Output(s):
%   sample_index    - [num array] sample indices corresponding to sample
%                     time
%   sample_yn       - [logical array] true: this sample index corresponding
%                     to physically collected data
%
% Note:
%
% See also SampleIndex2Time.

% Copyright 2019-2023 Richard J. Cui. Created: Sun 05/05/2019 10:29:21.071 PM
% $Revision: 1.5 $  $Date: Mon 10/09/2023 12:22:37.852 AM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% TODO: need more accurate estimate from uUTC to sample indexes
% parse inputs
% ------------
arguments
    this (1, 1) MultiscaleElectrophysiologyData
    sample_time (1, :) double
end % positional

arguments
    options.st_unit (1, 1) string {mustBeMember(options.st_unit, ["index", "uutc", "second", "minute", "hour", "day"])} = "index"
end %

st_unit = options.st_unit;

% TODO: convert to uUTC relative to the first sample of the recording
switch lower(st_unit)
    case 'msec'
        sample_time = round(sample_time * 1e3);
    case 'second'
        sample_time = round(sample_time * 1e6);
    case 'minute'
        sample_time = round(sample_time * 60 * 1e6);
    case 'hour'
        sample_time = round(sample_time * 60 * 60 * 1e6);
    case 'day'
        sample_time = round(sample_time * 24 * 60 * 60 * 1e6);
end % switch

% set paras
% ----------
sample_index = zeros(size(sample_time));
sample_yn = false(size(sample_time));
[sorted_st, orig_index] = sort(sample_time);

if isempty(this.ContinuityCorrected)
    cont = this.analyzeContinuity;
else
    cont = this.ContinuityCorrected;
end % if

% within continuous segment
% -------------------------
fs = this.ChannelMetadata.metadata.sampling_frequency;

% get relative start and end time in uUTC
x = cont{:, {'start_time', 'end_time'}}; % in uUTC
cont_start_end = x - x(1, 1); % in relative uUTC

% choose continuity segment that in the range of sample indexes
num_st = numel(sorted_st);
% TODO: need to consider the case when the start and end time are in the discontinuity
sel_cont_ind = sorted_st(1) <= cont_start_end(:, 2) ...
    & sorted_st(num_st) >= cont_start_end(:, 1);
sel_cont = cont(sel_cont_ind, :); % select the segment of continuity in the
% range of sorted_si
sel_cont_start_end = cont_start_end(sel_cont_ind, :);
[sorted_sample_index, sorted_sample_yn] = inContLoopCont(fs, sel_cont_start_end, ...
    sel_cont, sorted_st); % TODO

% output
% ------
sample_index(orig_index) = sorted_sample_index;
sample_yn(orig_index) = sorted_sample_yn;

end

% =========================================================================
% subroutines
% =========================================================================
function [s_index, s_yn] = inContLoopCont(fs, cont_se, cont, sorted_st)
% map sample time to sample index within continuous segment

arguments
    fs (1, 1) double % sampling frequency
    cont_se (:, 2) double % start and end time in uUTC
    cont (:, :) table % continuity table
    sorted_st (:, 1) double % sorted sample time in uUTC
end % positional

% check sample time one by one
% ----------------------------
num_st = numel(sorted_st); % number of sample time
s_index = nan(num_st, 1); % sample index
s_yn = false(num_st, 1); % sample yn

for k = 1:num_st
    st_k = sorted_st(k);
    
    % * check if the sample time is in the continuity segment
    idx_k = find(st_k >= cont_se(:, 1) & st_k <= cont_se(:, 2), 1, 'first');
    
    if isempty(idx_k)
        continue
    end % if
    
    s_yn(k) = true;
    
    % * get the sample index
    s_index(k) = round((st_k - cont_se(idx_k, 1)) .* fs / 1e6 ...
        + double(cont.start_index(idx_k)));
end % for

end % function

% [EOF]
