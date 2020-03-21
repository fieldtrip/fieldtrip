function bid = readBlockIndexData(this, varargin)
% MultiscaleElectrophysiologyFile_3p0.READBLOCKINDEXDATA read block index data from MEF 3.0 file
% 
% Syntax:
%   bid = readBlockIndexData(this)
%   bid = readBlockIndexData(__, channel)
% 
% Imput(s):
%   this            - [obj] MultiscaleElectrophysiologyFile_3p0 object
%   channel         - [struct] (opt) channel metadata structure
% 
% Output(s):
%   bid             - [N x 8 table] N is the number of blocks indexed. Each
%                     row has seven varialbes:
%                     Segment       : [num] the INDEX of segment of the 
%                                     recorded data block
%                     Block         : [num] the INDEX of the block in the 
%                                     segment
%                     FileOffset    : [num] offset in bytes of the 1st data
%                                     in the block in MEF 3.0 file
%                     StartTime     : [num] lower bound of the first sample
%                                     time recorded in the block (in uUTC)
%                     EndTime       : [num] upper bound of the last sample
%                                     time recorded in the block (in uUTC)
%                     StartSample   : [num] INDEX of the 1st sample of the
%                                     block in MEF file (1st sample index
%                                     in the file is zero in MEF file;
%                                     change to Matlab convention to start
%                                     at one)
%                     NumOfSamples  : [num] number of samples in the block
%                     REDBlockFlags : [Uint8] RED Block flags
%                                     Bit 0     - Discontinuity bit (1 dis.)
%                                     Bit 1     - Level 1 encrypted block
%                                                 bit
%                                     Bit 2     - Level 2 encrypted block
%                                                 bit
%                                     Bits 3-7  - reserved for future use
% 
% Note:
%   All indexes start at one (1), using matlab convention.
% 
%   See the details of MEF file at https://msel.mayo.edu/codes.html.
% 
% See also .

% Copyright 2020 Richard J. Cui. Created: Wed 02/05/2020 10:19:17.599 AM
% $Revision: 0.3 $  $Date: Wed 03/11/2020  4:41:01.147 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});
channel = q.channel;
if isempty(channel)
    channel = this.Channel;
end % if

% =========================================================================
% main
% =========================================================================
fs = this.ChanSamplingFreq;
MPS = this.MPS;
var_names = {'Segment', 'Block', 'FileOffset', 'StartTime', 'EndTime',...
    'StartSample', 'NumOfSamples', 'REDBlockFlags'};
bid = read_mef_bid(channel, MPS, fs, var_names);

% update
this.BlockIndexData = bid;
end

% =========================================================================
% subroutines
% =========================================================================
function bid = read_mef_bid(channel, MPS, fs, var_names)

num_seg = channel.number_of_segments;
% get number of blocks in each segments
num_blk = zeros(1, num_seg);
segments = channel.segments;
for k = 1:num_seg
    tsi_k = segments(k);
    num_blk(k) = numel(tsi_k.time_series_indices);    
end % for

% get the matrix
bid_mat = zeros(sum(num_blk), numel(var_names));
row_index = 0;
for j = 1:num_seg
    num_blk_j = num_blk(j);
    seg_j = segments(j);
    for k = 1:num_blk_j
        row_index = row_index+1;
        
        tsi_jk = seg_j.time_series_indices(k); % TSI of block k of segement j
        % segment index
        bid_mat(row_index, 1) = j;
        % block index
        bid_mat(row_index, 2) = k;
        % file_offset
        bid_mat(row_index, 3) = double(tsi_jk.file_offset);
        % start_time
        st_jk = double(tsi_jk.start_time);
        bid_mat(row_index, 4) = st_jk;
        % end_time
        num_samples_jk = double(tsi_jk.number_of_samples);
        et_jk = st_jk+num_samples_jk*MPS/fs;
        bid_mat(row_index, 5) = et_jk; % end_time
        % start_sample
        bid_mat(row_index, 6) = tsi_jk.start_sample+1; % change to matlab convention
        % number_of_samples
        bid_mat(row_index, 7) = num_samples_jk;
        % RED_block_flags
        bid_mat(row_index, 8) = tsi_jk.RED_block_flags;
    end % for
end % for

% get the BID table
bid = array2table(bid_mat, 'VariableNames', var_names);

end % function

function q = parseInputs(varargin)

% default
default_ch = struct([]);

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addOptional('channel', default_ch, @isstruct);

% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % funciton

% [EOF]