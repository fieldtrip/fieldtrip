function S = LoadSpikes(tfilelist,varargin)
%
% S = LoadSpikes(tfilelist)
%
% inp: tfilelist is a cellarray of strings, each of which is a
%   tfile to open.  Note: this is incompatible with version unix3.1.
% out: Returns a cell array such that each cell contains a ts 
%   object (timestamps which correspond to times at which the cell fired)

% ADR 1998
%  version L4.0
%  status: PROMOTED

tsflag = 'sec';
display = 1;

try
extract_varargin;
end

%-------------------
% Check input type
%-------------------
if ~isa(tfilelist, 'cell')
   error('LoadSpikes: tfilelist should be a cell array.');
end

nFiles = length(tfilelist);

%--------------------
% Read files
%--------------------

%fprintf(2, 'Reading %d files.', nFiles);

% for each tfile
% first read the header, the read a tfile 
% note: uses the bigendian modifier to ensure correct read format.

S = cell(nFiles, 1);
for iF = 1:nFiles
%   if display
%       DisplayProgress(iF, nFiles, 'Title', 'LoadSpikes');
%   end
    tfn = tfilelist{iF};
    if ~isempty(tfn)
        tfp = fopen(tfn, 'rb','ieee-be');
        if (tfp == -1)
            warning([ 'Could not open tfile ' tfn]);
        end
        
        ReadHeader(tfp);    
        S{iF} = fread(tfp,inf,'uint32');    %read as 32 bit ints
        
        % Set appropriate time units.
        switch (tsflag)
            case 'sec'
                S{iF} = S{iF}/10000;
            case 'ts'
                S{iF} = S{iF};
            case 'ms'
                S{iF} = S{iF}*10000*1000;
            otherwise
                error('LoadSpikes called with invalid tsflag.');
        end
        
%       S{iF} = ts(S{iF});
        
        fclose(tfp);
        
    end         % if tfn valid
end     % for all files
fprintf(2,'\n');
