function seg_cont = analyzeContinuity(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE_2P1.ANALYZECONTINUITY Analyze continuity of sampling in MEF 2.1 data
% 
% Syntax:
%   seg_cont = analyzeContinuity(this)
%   seg_cont = analyzeContinuity(this, wholename)
%   seg_cont = analyzeContinuity(this, wholename, password)
% 
% Imput(s):
%   this            - [obj] MultiscaleElectrophysiologyFile object
%   wholename       - [str] filepath + filename of MEF file
%   password        - [str] subject password of the data
% 
% Output(s):
%   seg_cont        - [table] N x 7, information of segments of continuity
%                     of sampling in the data file.  The 7 variable names
%                     are:
%                     BlockStart        : [num] index of start data block
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
%                                         in this segment of continuity
%                                         sampling
% 
% Note:
%   See the details of MEF file at https://github.com/benbrinkmann/mef_lib_2_1
% 
% See also .

% Copyright 2019-2020 Richard J. Cui. Created: Sat 05/04/2019 10:35:40.540 PM
% $Revision: 0.4 $  $Date: Sun 02/09/2020 10:40:13.746 AM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, varargin{:});

if isempty(q)
    wholename = fullfile(this.FilePath, this.FileName);
    pw = this.SubjectPassword;
else
    wholename = fullfile(q.filepath, q.filename);
    pw = q.password;
    % update
    this.FilePath = q.filepath;
    this.FileName = q.filename;
    this.Password = q.password;
end % if

% get header of file
% ------------------
if isempty(this.Header)
    this.readHeader(wholename, pw);
end % if
header = this.Header;

% table of segment of continuity
% ------------------------------
% number of continuity segment = number of indexed discontinuity blocks
num_seg_cont = header.number_of_discontinuity_entries; 
varNames = {'BlockStart', 'BlockEnd', 'SampleTimeStart', 'SampleTimeEnd',...
            'SampleIndexStart', 'SampleIndexEnd', 'SegmentLength'};
seg_cont = table('Size', [num_seg_cont, numel(varNames)], 'VariableTypes',...
    {'double', 'double', 'double', 'double', 'double', 'double', 'double'},...
    'VariableNames', varNames);

% read discontinuity indices
fp = fopen(wholename, 'r');
if fp < 0, return; end % if

fseek(fp, header.discontinuity_data_offset, 'bof');
a = fread(fp, header.number_of_discontinuity_entries, 'uint64');
a = a(:);
fclose(fp);

% block start-end
if a == 0
    blk_start = 1;
    blk_end = header.number_of_index_entries;
else
    blk_start = a; % 1st block indexed as one
    blk_end = [a(2:end)-1; header.number_of_index_entries];
end % if
seg_cont{:, {'BlockStart', 'BlockEnd'}} = [blk_start, blk_end];

% other info
if isempty(this.BlockIndexData)
    this.readBlockIndexData;
end % if
bid = this.BlockIndexData;
fs = header.sampling_frequency;
MPS = 1e6; % microseconds per second
wh = waitbar(0,'Analyzing signal continuity...');
fprintf('Analyzing signal continuity...')
for k = 1:num_seg_cont
    waitbar(k/num_seg_cont)
    
    blk_start_k = seg_cont.BlockStart(k);
    % sample time start
    seg_cont.SampleTimeStart(k) = bid.SampleTime(blk_start_k);
    % sample index start
    seg_cont.SampleIndexStart(k) = bid.SampleIndex(blk_start_k);
    
    blk_end_k = seg_cont.BlockEnd(k);
    % sample time end
    end_blk_header_k = this.readBlockHeader(blk_end_k);
    end_blk_start_time = bid.SampleTime(blk_end_k);
    blk_len_k = end_blk_header_k.block_length; % total number of samples
    seg_cont.SampleTimeEnd(k) = end_blk_start_time+round(MPS/fs)*blk_len_k-1;
    % sample index end
    end_blk_start_ind = bid.SampleIndex(blk_end_k);
    seg_cont.SampleIndexEnd(k) = end_blk_start_ind+blk_len_k-1;
    
    % segment length
    seg_cont.SegmentLength(k) = seg_cont.SampleIndexEnd(k)-seg_cont.SampleIndexStart(k)+1;
end % for
fprintf('Done!\n')
close(wh)

this.Continuity = seg_cont;
end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(this, varargin)

% defaults
default_pw = this.SubjectPassword;

% parse rules
p = inputParser;
p.addOptional('wholename', '', @isstr);
p.addOptional('password', default_pw, @isstr);

% parse and return the results
p.parse(varargin{:});
if isempty(p.Results.wholename)
    q = [];
else
    [fp, fn, ext] = fileparts(p.Results.wholename);
    q.filepath = fp;
    q.filename = [fn, ext];
    q.password = p.Results.password;
end % if

end % function

% [EOF]

