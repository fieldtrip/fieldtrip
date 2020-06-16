function bid = readBlockIndexData(this, varargin)
% MULTISCALEELECTROPHYSIOLOGYFILE_2P1.READBLOCKINDEXDATA Read Block Index Data from MEF 2.1 file
% 
% Syntax:
%   bid = readBlockIndexData(this)
%   bid = readBlockIndexData(this, wholename)
%   bid = readBlockIndexData(this, wholename, password)
% 
% Imput(s):
%   this            - [obj] MultiscaleElectrophysiologyFile_2p1 object
%   wholename       - [str] filepath + filename of MEF file
%   password        - [str] password of the data
% 
% Output(s):
%   bid             - [N x 3 table] N is the number of blocks indexed. Each
%                     row has three varialbes (triplets):
%                     SampleTime    : [ui8] first sample time recorded in
%                                     the block (in uUTC)
%                     FileOffset    : [ui8] offset in bytes of the 1st data
%                                     in the block in MEF file (including
%                                     block header bytes)
%                     SampleIndex   : [ui8] index of the 1st sample of the
%                                     block in MEF file (1st sample index
%                                     in the file is zero in MEF file;
%                                     change to Matlab convention to start
%                                     with one)
% 
% Note:
%   See the details of MEF file at https://github.com/benbrinkmann/mef_lib_2_1
% 
% See also .

% Copyright 2019-2020 Richard J. Cui. Created: Tue 04/30/2019 10:11:41.380 PM
% $Revision: 0.3 $  $Date: Wed 02/05/2020 10:38:49.323 AM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

q = parseInputs(this, varargin{:});

if ~isempty(q)
    % update
    this.FilePath = q.filepath;
    this.FileName = q.filename;
    this.Password = q.password;
    this.readHeader;
end % if

varNames = {'SampleTime', 'FileOffset', 'SampleIndex'};
bid = read_mef_bid(this, varNames);
bid.SampleIndex = bid.SampleIndex+1; % change to Matlab convention

this.BlockIndexData = bid;

end

% =========================================================================
% subroutines
% =========================================================================
function bid = read_mef_bid(this, varNames)
% read block indices data from MEF

% check header of MEF
% -------------------
if isempty(this.Header)
    this.readHeader;
end % if
header = this.Header;

% read bid
% --------
wholename = fullfile(this.FilePath, this.FileName);
fp = fopen(wholename, 'r');
if fp < 0, return; end

fseek(fp, header.index_data_offset, 'bof');
a = fread(fp, [3, header.number_of_index_entries], 'uint64', 0, 'l');
bid = array2table(a.', 'VariableNames', varNames);

fclose(fp);

end % function

function q = parseInputs(this, varargin)

% defaults
default_pw = this.SessionPassword; % use session pw

% parse rules
p = inputParser;
p.addOptional('wholename', '', @isstr);
p.addOptional('password', default_pw, @isstr);

% parse and return the results
p.parse(varargin{:});
wholename = p.Results.wholename;
if isempty(wholename)
    q = [];
else
    if exist(wholename, 'file') ~= 2
        error('Cannot find %s.', wholename);
    end % if

    [fp, fn, ext] = fileparts(wholename);
    q.filepath = fp;
    q.filename = [fn, ext];
    q.password = p.Results.password;
end % if

end % function

% [EOF]

