function [sample_time, sample_yn] = SampleIndex2Time(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE.SAMPLEINDEX2TIME Convert sample index to sample time
% 
% Syntax
%   [sample_time, sample_yn] = SampleIndex2Time(this, sample_index)
%   [sample_time, sample_yn] = SampleIndex2Time(__, st_unit)
% 
% Input(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
%   sample_index    - [num array] array of sample index (must be integers)
%   st_unit         - [str] (optional) sample time unit: 'uUTC' (default)
%                     or 'u', 'mSec', 'Second' or 's', 'Minute' or 'm', 'Hour' or
%                     'h' and 'Day' or 'd'.
% 
% Output(s):
%   sample_time     - [num array] sample time corresponding to sample
%                     indices (default unit: uUTC)
%   sample_yn       - [logical array] true: this sample time corresponding
%                     to physically collected data
% 
% Note:
%   An error less than one sample time may occure.
% 
% See also SampleTime2Index.

% Copyright 2019-2020 Richard J. Cui. Created: Mon 05/06/2019  9:29:08.940 PM
% $Revision: 1.1 $  $Date: Wed 02/05/2020  9:10:33.935 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% parse inputs
% ------------
q = parseInputs(this, varargin{:});
sample_index = q.sample_index;
st_unit = q.st_unit;

% set paras
% ----------
sample_time = zeros(size(sample_index));
sample_yn = false(size(sample_index));
[sorted_si, orig_index] = sort(sample_index);

if isempty(this.Continuity)
    this.analyzeContinuity;
end % if
cont = this.Continuity;
number_of_discontinuity_entries = height(cont);

% within continuous segment
% -------------------------
cont_start_end = cont{:, {'SampleIndexStart', 'SampleIndexEnd'}};

% choose continuity segment that in the range of sample indexes
num_si = numel(sorted_si);
sel_cont_ind = sorted_si(1) <= cont_start_end(:, 2) ...
    & sorted_si(num_si) >= cont_start_end(:, 1);
sel_cont = cont(sel_cont_ind, :); % select the segment of continuity in the
                                  % range of sorted_si
sel_cont_start_end = cont_start_end(sel_cont_ind, :);                                  
[sorted_sample_time, sorted_sample_yn] = inContLoopCont(sel_cont_start_end,...
    sel_cont, sorted_si);

% within discontinous segment
% ----------------------------
if number_of_discontinuity_entries > 1 % if discont in recording
    a = cont_start_end.';
    b = cat(1, -inf, a(:), inf);
    discont_start_end = reshape(b, 2, numel(b)/2).';
    not_discont_index = discont_start_end(:, 2)-discont_start_end(:, 1) == 1;
    discont_start_end(not_discont_index, :) = [];
    discont_index = ~not_discont_index;
    sel_cont_ind = [true; discont_index(2:end-2); true];
    sel_cont = cont(sel_cont_ind, :);
    [sorted_sample_time, sorted_sample_yn] = inDiscontLoopDiscont(...
        sorted_sample_time, sorted_sample_yn, discont_start_end, sel_cont,...
        sorted_si);
end % if

% output
% ------
switch lower(st_unit)
    case 'msec'
        sorted_sample_time = sorted_sample_time/1e3;
    case 'second'
        sorted_sample_time = sorted_sample_time/1e6;
    case 'minute'
        sorted_sample_time = sorted_sample_time/(60*1e6);
    case 'hour'
        sorted_sample_time = sorted_sample_time/(60*60*1e6);
    case 'day'
        sorted_sample_time = sorted_sample_time/(24*60*60*1e6);
end % switch
sample_time(orig_index) = sorted_sample_time;
sample_yn(orig_index) = sorted_sample_yn;

end

% =========================================================================
% subroutines
% =========================================================================
function [s_time, s_yn] = inContLoopCont(cont_se, cont, sorted_si)

blk_len = 2^20; % length of one block of sample indexed
num_si = numel(sorted_si);
if num_si >= 3*blk_len % set verbose
    verb = true;
else
    verb = false;
end % if

s_time = zeros(size(sorted_si));
s_yn = false(size(sorted_si));

num_blk = ceil(num_si/blk_len);
if verb, wh = waitbar(0, 'Coverting sample indexes to time...'); end % if
for k = 1:num_blk
    start_k = (k-1)*blk_len+1;
    end_k = k*blk_len;
    if end_k >= num_si
        end_k = num_si;
    end % if
    si_k = sorted_si(start_k:end_k);
    [s_time_k, s_yn_k] = inContLoopCont_blk(cont_se, cont, si_k);
    
    s_time(start_k:end_k) = s_time_k;
    s_yn(start_k:end_k) = s_yn_k;
    if verb, waitbar(k/num_blk, wh), end % if
end % for
if verb, close(wh), end % if

end % function

function [s_time, s_yn] = inDiscontLoopDiscont(s_time, s_yn, discont_se,...
    cont, sorted_si)

MPS = 1e6;
num_seg = size(discont_se, 1); % number of segments
for k = 1:num_seg
    start_k = discont_se(k, 1);
    end_k = discont_se(k, 2);
    ind_k = sorted_si > start_k & sorted_si < end_k;
    
    if sum(ind_k) ~= 0
        if start_k == -inf
            st_k = cont.SampleTimeStart(k);
            index_diff = sorted_si(ind_k)-end_k;
            time_diff = index_diff*MPS/fs;
            sorted_st_k = st_k+time_diff;
        else
            st_k = cont.SampleTimeEnd(k-1);
            index_diff = sorted_si(ind_k) - start_k;
            time_diff = (index_diff-1)*MPS/fs;
            sorted_st_k = st_k+time_diff+1;
        end % if
        s_time(ind_k) = round(sorted_st_k); % uUTC integer
        s_yn(ind_k) = false;
    end % if
end % for
end % function

function [s_time, s_yn] = inContLoopCont_blk(cont_se, cont, sorted_si)
% within continuous segment loop through continuity segments

s_time = zeros(size(sorted_si));
s_yn = false(size(sorted_si));

num_seg = size(cont_se, 1); % number of segments
for k = 1:num_seg
    start_k = cont_se(k, 1);
    end_k = cont_se(k, 2);
    ind_k = sorted_si >= start_k & sorted_si <= end_k;
    
    if sum(ind_k) ~= 0
        st_k = cont.SampleTimeStart(k);
        et_k = cont.SampleTimeEnd(k);
        slop_k = (et_k-st_k)/(end_k-start_k);
        index_diff = sorted_si(ind_k)-start_k;
        sorted_st_k = st_k+slop_k*index_diff; % use linear approximation
        s_time(ind_k) = round(sorted_st_k); % uUTC integer
        s_yn(ind_k) = true;
    end % if
end % for

end % function

function q = parseInputs(varargin)

% defaults
defaultSTUnit = 'uutc';
expectedSTUnit = {'uutc', 'msec', 'second', 'minute', 'hour', 'day'};

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('sample_index', @isnumeric);
p.addOptional('st_unit', defaultSTUnit,...
    @(x) any(validatestring(x, expectedSTUnit)));

% parse and return the results
p.parse(varargin{:});
q.this = p.Results.this;
q.sample_index = p.Results.sample_index;
q.st_unit = p.Results.st_unit;

end % function

% [EOF]