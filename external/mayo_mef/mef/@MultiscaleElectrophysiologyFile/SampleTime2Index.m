function [sample_index, sample_yn] = SampleTime2Index(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE.SAMPLETIME2INDEX Convert sample time to sample index
% 
% Syntax:
%   [sample_index, sample_yn] = SampleTime2Index(this, sample_time)
%   [sample_index, sample_yn] = SampleTime2Index(__, st_unit)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
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

% Copyright 2019-2020 Richard J. Cui. Created: Sun 05/05/2019 10:29:21.071 PM
% $Revision: 0.7 $  $Date: Wed 02/05/2020  8:39:11.992 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% TODO: need more accurate estimate from uUTC to sample indexes
% parse inputs
% ------------
q = parseInputs(this, varargin{:});
sample_time = q.sample_time;
st_unit = q.st_unit;
switch lower(st_unit) % convert to uUTC
    case 'msec'
        sample_time = round(sample_time*1e3);
    case 'second'
        sample_time = round(sample_time*1e6);
    case 'minute'
        sample_time = round(sample_time*60*1e6);
    case 'hour'
        sample_time = round(sample_time*60*60*1e6);
    case 'day'
        sample_time = round(sample_time*24*60*60*1e6);
end % switch

% set paras
% ----------
sample_index = zeros(size(sample_time));
sample_yn = false(size(sample_time));
[sorted_st, orig_index] = sort(sample_time);

if isempty(this.Continuity)
    cont = this.analyzeContinuity;
else
    cont = this.Continuity;
end % if
number_of_discontinuity_entries = height(cont);

% within continuous segment
% -------------------------
cont_start_end = cont{:, {'SampleTimeStart', 'SampleTimeEnd'}};

% choose continuity segment that in the range of sample indexes
num_st = numel(sorted_st);
sel_cont_ind = sorted_st(1) <= cont_start_end(:, 2)...
    & sorted_st(num_st) >= cont_start_end(:, 1);
sel_cont = cont(sel_cont_ind, :); % select the segment of continuity in the
                                  % range of sorted_si
sel_cont_start_end = cont_start_end(sel_cont_ind, :);                                  
[sorted_sample_index, sorted_sample_yn] = inContLoopCont(sel_cont_start_end,...
    sel_cont, sorted_st);

% within discontinous segment
% ----------------------------
% TODO: need to segment the data for very large chunck of times in the
% discontinuity segments, similar with the case of continuity
if number_of_discontinuity_entries > 1 % if discont in recording
    a = cont_start_end.';
    b = cat(1, -inf, a(:), inf);
    discont_start_end = reshape(b, 2, numel(b)/2).';
    num_seg = size(discont_start_end, 1); % number of segments
    for k = 1:num_seg
        start_k = discont_start_end(k, 1);
        end_k = discont_start_end(k, 2);
        ind_k = sorted_st > start_k & sorted_st < end_k;
        
        if sum(ind_k) ~= 0
            sorted_sample_index(ind_k) = NaN; % if between one index difference
            sorted_sample_yn(ind_k) = false;
        end % if
    end % for
end % if

% output
% ------
sample_index(orig_index) = sorted_sample_index;
sample_yn(orig_index) = sorted_sample_yn;

end

% =========================================================================
% subroutines
% =========================================================================
function [s_index, s_yn] = inContLoopCont(cont_se, cont, sorted_st)

blk_len = 2^20; % length of one block of sample times
num_st = numel(sorted_st);
if num_st >= 3*blk_len
    verb = true;
else
    verb = false;
end % if

s_index = zeros(size(sorted_st));
s_yn = false(size(sorted_st));

num_blk = ceil(num_st/blk_len);
if verb, wh = waitbar(0, 'Coverting sample times to indexes...'); end % if
for k = 1:num_blk
    start_k = (k-1)*blk_len+1;
    end_k = k*blk_len;
    if end_k >= num_st
        end_k = num_st;
    end % if
    st_k = sorted_st(start_k:end_k);
    [s_index_k, s_yn_k] = inContLoopCont_blk(cont_se, cont, st_k);
    
    s_index(start_k:end_k) = s_index_k;
    s_yn(start_k:end_k) = s_yn_k;
    if verb, waitbar(k/num_blk, wh), end % if
end % for
if verb, close(wh), end % if

end % function

function [s_index, s_yn] = inContLoopCont_blk(cont_se, cont, sorted_st)
% within continuous segment loop through continuity segments

s_index = zeros(size(sorted_st));
s_yn = false(size(sorted_st));

num_seg = size(cont_se, 1); % number of segments
for k = 1:num_seg
    start_k = cont_se(k, 1); % start time
    end_k = cont_se(k, 2); % end time
    ind_k = sorted_st >= start_k & sorted_st <= end_k;
    
    if sum(ind_k) ~= 0
        si_k = cont.SampleIndexStart(k);
        ei_k = cont.SampleIndexEnd(k);
        slop_k = (ei_k-si_k)/(end_k-start_k);
        time_diff = sorted_st(ind_k) - start_k;
        sorted_ti_k = si_k+slop_k*time_diff;
        s_index(ind_k) = round(sorted_ti_k); % align time to the nearest index
        s_yn(ind_k) = true;
    end % if
end % for

end % funciton

function q = parseInputs(varargin)

% defaults
defaultSTUnit = 'uutc';
expectedSTUnit = {'uutc', 'msec', 'second', 'minute', 'hour', 'day'};

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('sample_time', @isnumeric);
p.addOptional('st_unit', defaultSTUnit,...
    @(x) any(validatestring(x, expectedSTUnit)));

% parse and return the results
p.parse(varargin{:});
q.this = p.Results.this;
q.sample_time = p.Results.sample_time;
q.st_unit = p.Results.st_unit;

end % function

% [EOF]
