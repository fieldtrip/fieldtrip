function blk_header = readBlockHeader(this, BlockIndex)
% READBLOCKHEADER Read header of the specified data block
% 
% Syntax:
%   blk_header = readBlockHeader(this, BlockIndex)
% 
% Imput(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
%   BlockIndex      - [1 x N array] N block indices (1st block indexed as
%                     one)
% 
% Output(s):
%   blk_header      - [1 x N struct array] contains N block header
%                     structure:
%                     .crc      : [ui4] Cyclically Redundant Checksum (of
%                                 the data block)
%                     .compressed_block_length
%                               : [ui4] number of bytes in the compressed
%                                 block, excluding block header
%                     .block_start_time
%                               : [ui8] start recording time (uUTC time)
%                     .difference_length
%                               : [ui4] difference data length in bytes
%                     .block_length
%                               : [ui4] number of data samples encoded in
%                                 the block
%                     .maximum_data_value
%                               : [si3] the maximum raw value encoded in
%                                 the data block
%                     .minimum_data_value
%                               : [si3] the minimum raw value encoded in
%                                 the data block
%                     .block_flag
%                               : [ui1] Bit 0: 0 no discontinuity; 1 this
%                                 block began after a discontinuity in
%                                 recording, or the first block in a file.
%                                 Bits 1-7 unused.
%                     .block_statistics
%                               : [ui1] statistical model of difference
%                                 value for the block
% 
% Note:
%   See the details of MEF file at https://github.com/benbrinkmann/mef_lib_2_1
% 
% See also .

% Copyright 2019 Richard J. Cui. Created: Fri 05/03/2019  6:56:26.566 PM
% $Revision: 0.1 $  $Date: Fri 05/03/2019  6:56:26.566 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

q = parseInputs(this, BlockIndex);
blk_index = q.BlockIndex;

wholename = fullfile(this.FilePath, this.FileName);
if exist(wholename, 'file') ~= 2
    error('Cannot find file %s', wholename)
end % if

if isempty(this.BlockIndexData)
    this.readBlockIndexData;
end % if
bid = this.BlockIndexData;

fp = fopen(wholename, 'r');
if fp < 0, return; end % if

blk_header = struct([]);
for bi_k = blk_index
    fseek(fp, bid.FileOffset(bi_k), 'bof'); % move the pointer to block start
    blk_header_k = read_mef_blkheader(fp);
    blk_header = cat(2, blk_header, blk_header_k);
end % for

fclose(fp);

end

% =========================================================================
% subroutines
% =========================================================================
function r = read24bits(a)
% a         - 24 bits signed integer
% r         - int32

b = flipud(a);
c = dec2bin(b, 8)';
d = c(:)';
e = [d(ones(1, 8)), d]; % extend the sign
r = typecast(uint32(bin2dec(e)), 'int32');

end % function

function blkh = read_mef_blkheader(fp)
% read one block header
%
% fp        - file pointer to the beginning of the block

blkh.crc                    = fread(fp, 1, 'uint32');
blkh.compressed_block_length= fread(fp, 1, 'uint32');
blkh.block_start_time       = fread(fp, 1, 'uint64');
blkh.difference_length      = fread(fp, 1, 'uint32');
blkh.block_length           = fread(fp, 1, 'uint32');
blkh.maximum_data_value     = read24bits(fread(fp, 3, 'uchar')); % read 24 bits, not sure
blkh.minimum_data_value     = read24bits(fread(fp, 3, 'uchar')); % note sure
blkh.block_flag             = fread(fp, 1, 'uint8');
blkh.block_statistics       = fread(fp, 256, 'uint8');

end % function

function q = parseInputs(varargin)

% defaults

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('BlockIndex', @isnumeric);

% parse and return the results
p.parse(varargin{:});
q.this = p.Results.this;
q.BlockIndex = p.Results.BlockIndex;

end % function

% [EOF]

