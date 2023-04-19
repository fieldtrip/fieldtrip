function [ compdata ] = fiff_read_ctf_comp(fid,node,chs,ch_rename)

%
% [ compdata ] = fiff_read_ctf_comp(fid,node,chs,ch_rename)
%
% Read the CTF software compensation data from the given node
%


%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.8  2007/11/13 10:55:32  msh
%   Specify the second argument to all calls to the exist function
%
%   Revision 1.7  2006/09/25 19:48:16  msh
%   Added projection item kinds to fiff_define_constants
%   Changed some fields to int32 in surface structures
%
%   Revision 1.6  2006/08/09 15:22:51  msh
%   Removed debug printout.
%
%   Revision 1.5  2006/06/22 21:22:46  msh
%   Take into account the possibility of calibrated compensation matrices
%
%   Revision 1.4  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.3  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.2  2006/04/12 10:29:02  msh
%   Made evoked data writing compatible with the structures returned in reading.
%
%   Revision 1.1  2006/04/10 23:26:54  msh
%   Added fiff reading routines
%
%

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

me='MNE:fiff_read_ctf_comp';
if nargin == 3
    ch_rename = {};
elseif nargin ~= 4
    error(me,'Incorrect number of arguments');
end

compdata = struct('ctfkind',{},'kind',{},'save_calibrated',{},'rowcals',{},'colcals',{},'data',{});
comps = fiff_dir_tree_find(node,FIFF.FIFFB_MNE_CTF_COMP_DATA);

for k = 1:length(comps)
    node = comps(k);
    %
    %   Read the data we need
    %
    mat  = fiff_read_named_matrix(fid,node,FIFF.FIFF_MNE_CTF_COMP_DATA);
    for p = 1:node.nent
        kind = node.dir(p).kind;
        pos  = node.dir(p).pos;
        if kind == FIFF.FIFF_MNE_CTF_COMP_KIND
            tag = fiff_read_tag(fid,pos);
            break;
        end
    end
    if ~exist('tag','var')
        error(me,'Compensation type not found');
    end
    %
    %   Get the compensation kind and map it to a simple number
    %
    one.ctfkind = tag.data;
    clear('tag');
    one.kind    = int32(-1);
    if one.ctfkind     == hex2dec('47314252')
        one.kind = int32(1);
    elseif one.ctfkind == hex2dec('47324252')
        one.kind = int32(2);
    elseif one.ctfkind == hex2dec('47334252')
        one.kind = int32(3);
    else
        one.kind = one.ctfkind;
    end
    for p = 1:node.nent
        kind = node.dir(p).kind;
        pos  = node.dir(p).pos;
        if kind == FIFF.FIFF_MNE_CTF_COMP_CALIBRATED
            tag = fiff_read_tag(fid,pos);
            break;
        end
    end
    if ~exist('tag','var')
        calibrated = false;
    else
        calibrated = tag.data;
    end
    one.save_calibrated = calibrated;
    one.rowcals = ones(1,size(mat.data,1));
    one.colcals = ones(1,size(mat.data,2));
    if ~calibrated
        %
        %   Calibrate...
        %
        %
        %   Do the columns first
        %
        for p  = 1:length(chs)
            ch_names{p} = chs(p).ch_name;
        end
        for col = 1:size(mat.data,2)
            p = strmatch(mat.col_names{col},ch_names,'exact');
            if isempty(p)
                error(me,'Channel %s is not available in data',mat.col_names{col});
            elseif length(p) > 1
                error(me,'Ambiguous channel %s',mat.col_names{col});
            end
            col_cals(col) = 1.0/(chs(p).range*chs(p).cal);
        end
        %
        %    Then the rows
        %
        for row = 1:size(mat.data,1)
            p = strmatch(mat.row_names{row},ch_names,'exact');
            if isempty(p)
                error(me,'Channel %s is not available in data',mat.row_names{row});
            elseif length(p) > 1
                error(me,'Ambiguous channel %s',mat.row_names{row});
            end
            row_cals(row) = chs(p).range*chs(p).cal;
        end
        mat.data            = diag(row_cals)*mat.data*diag(col_cals);
        one.rowcals         = row_cals;
        one.colcals         = col_cals;
    end
    one.data       = mat;
    compdata(k)    = one;
    clear('row_cals');
    clear('col_cals');
    one = fiff_rename_comp(one, ch_rename);
end

if length(compdata) > 0
    fprintf(1,'\tRead %d compensation matrices\n',length(compdata));
end

return;

end
