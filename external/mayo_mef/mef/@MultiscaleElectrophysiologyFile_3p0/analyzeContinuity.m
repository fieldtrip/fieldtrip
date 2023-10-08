function seg_cont = analyzeContinuity(this, varargin)
% MultiscaleElectrophysiologyFile_3p0.analyzeContinuity analyze continuity of sampling in MEF 3.0 data
%   
% Syntax:
%   seg_cont = analyzeContinuity(this)
%   seg_cont = analyzeContinuity(__, bid)
% 
% Imput(s):
%   this            - [obj] MultiscaleElectrophysiologyFile_3p0 object
%   bid             - [table] (opt) block index data (see readBlockIndexData.m
%                     for the detail) (default = this.BlockIndexData)
% 
% Output(s):
%   seg_cont        - [table] N x 9, information of segments of continuity
%                     of sampling in the data file.  The 9 variable names
%                     are:
%                     SegmentStart      : [num] index of start data segment
%                     BlockStart        : [num] index of start data block 
%                     SegmentEnd        : [num] index of end data segment
%                     BlockEnd          : [num] index of end data block
%                     SampleTimeStart   : [num] sample time of start
%                                         recording (in uUTC)
%                     SampleTimeEnd     : [num] sample time of end
%                                         recording (in uUTC)
%                     SampleIndexStart  : [num] sample index of start
%                                         recording
%                     SampleIndexEnd    : [num] sample index of end
%                                         recording
%                     SegmentLength     : [num] the total number of samples
%                                         in this segment of continuous
%                                         sampling
% 
% Note:
%   See the details of MEF 3.0 file at https://msel.mayo.edu/codes.html
% 
% See also readBlockIndexData.

% Copyright 2020 Richard J. Cui. Created: Wed 02/05/2020 10:19:17.599 AM
% $Revision: 0.5 $  $Date: Wed 04/19/2023 12:22:52.916 AM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
bid = q.bid;
if isempty(bid) && isempty(this.BlockIndexData)
    bid = this.readBlockIndexData;
elseif isempty(bid) && ~isempty(this.BlockIndexData)
    bid = this.readBlockIndexData;
end % if

% =========================================================================
% main
% =========================================================================

% find out index of continous chunks
% -----------------------------------
flag_2nd_check = false;

% (1) use time criterion
num_blks = height(bid); % total number of blocks
st = bid.StartTime(2:num_blks);
et = bid.EndTime(1:num_blks-1);
x = find(abs(et-st) > 1, 1); x = x(:); % make sure to get vertical vector
if isempty(x)
    chunk_index = [1, num_blks];
else
    bg_ind = [1; x+1]; % begin index
    ed_ind = [x; num_blks]; % end index
    chunk_index = [bg_ind, ed_ind];
end % if

if flag_2nd_check == true
    % (2) use RED_block_flags
    rbf = bid.REDBlockFlags; % RED block flags
    cont_ind = bitget(rbf, 1) == 1; % continuity index
    bg_ind_rbf = find(cont_ind); % begin index
    end_ind_rbf = find([cont_ind(2:end); true]); % end index
    chunk_index_rbf = [bg_ind_rbf, end_ind_rbf];
    
    % check consistency
    if sum(chunk_index-chunk_index_rbf) ~= 0
        warning('off', 'backtrace')
        warning('MultiscaleElectrophysiologyFile_3p0:analyzeContinuity:noConsistent',...
            'Continuity blocks are not consistent in continuity checking')
    end % if
end % if

% get the continuity table
% ------------------------
var_names = {'SegmentStart', 'BlockStart', 'SegmentEnd', 'BlockEnd', 'SampleTimeStart',...
    'SampleTimeEnd', 'SampleIndexStart', 'SampleIndexEnd', 'SegmentLength'};
num_chunk_cont = size(chunk_index, 1); % number of continuity chunks
chunk_cont = zeros(num_chunk_cont, numel(var_names));
for k = 1:num_chunk_cont
    bg_k = chunk_index(1);
    ed_k = chunk_index(2);
    
    % SegmentStart
    chunk_cont(k, 1) = bid.Segment(bg_k);
    % BlockStart
    chunk_cont(k, 2) = bid.Block(bg_k);
    % SegmentEnd
    chunk_cont(k, 3) = bid.Segment(ed_k);
    % BlockEnd
    chunk_cont(k, 4) = bid.Block(ed_k);
    % SampleTimeStart
    chunk_cont(k, 5) = bid.StartTime(bg_k);
    % SampleTimeEnd
    chunk_cont(k, 6) = bid.EndTime(ed_k);
    % SampleIndexStart
    chunk_cont(k, 7) = bid.StartSample(bg_k);
    % SampleIndexEnd
    seg_len = sum(bid.NumOfSamples(bg_k:ed_k)); % segment length
    chunk_cont(k, 8) = chunk_cont(k, 7)+seg_len-1;
    % SegmentLength
    chunk_cont(k, 9) = seg_len;
end % for
seg_cont = array2table(chunk_cont, 'VariableNames', var_names);

% update
% -------
this.Continuity = seg_cont;

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults
default_bid = table;

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('bid', default_bid, @istable);

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]