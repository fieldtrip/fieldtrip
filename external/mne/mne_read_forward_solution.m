function [fwd] = mne_read_forward_solution(fname,force_fixed,surf_ori,include,exclude)
%
% [fwd] = mne_read_forward_solution(fname,force_fixed,surf_ori,include,exclude)
%
% A forward solution from a fif file
%
% fname        - The name of the file
% force_fixed  - Force fixed source orientation mode? (optional)
% surf_ori     - Use surface based source coordinate system? (optional)
% include      - Include these channels (optional)
% exclude      - Exclude these channels (optional)
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.11  2008/06/07 21:22:10  msh
%   Added FIFF_FORWARD_SOLUTION_GRAD constant
%
%   Revision 1.10  2008/06/07 21:19:00  msh
%   Added reading of the field gradients.
%
%   Revision 1.9  2006/11/30 05:43:29  msh
%   Fixed help text in fiff_dir_tree_find
%   Fixed check for the existence of parent MRI block in mne_read_forward_solution
%
%   Revision 1.8  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%   Revision 1.7  2006/05/03 19:09:03  msh
%   Fixed even more compatibility issues.
%
%   Revision 1.6  2006/04/24 00:07:24  msh
%   Fixed error in the computation of the surface-based coordinate system.
%
%   Revision 1.5  2006/04/23 15:29:41  msh
%   Added MGH to the copyright
%
%   Revision 1.4  2006/04/21 21:33:07  msh
%   Fixed error in combination of the forward solution.
%   Fixed help text in mne_read_inverse_operator
%
%   Revision 1.3  2006/04/20 21:49:38  msh
%   Added mne_read_inverse_operator
%   Changed some of the routines accordingly for more flexibility.
%
%   Revision 1.2  2006/04/18 23:21:22  msh
%   Added mne_transform_source_space_to and mne_transpose_named_matrix
%
%   Revision 1.1  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%

me='MNE:mne_read_forward_solution';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin == 1
    include = [];
    exclude = [];
    force_fixed = false;
    surf_ori = false;
elseif nargin == 2
    include = [];
    exclude = [];
    surf_ori = false;
elseif nargin == 3
    include = [];
    exclude = [];
elseif nargin == 4
    exclude = [];
elseif nargin ~= 5
    error(me,'Incorrect number of arguments');
end
%
%   Open the file, create directory
%
fprintf(1,'Reading forward solution from %s...\n',fname);
[ fid, tree ] = fiff_open(fname);
%
%   Find all forward solutions
%
fwds = fiff_dir_tree_find(tree,FIFF.FIFFB_MNE_FORWARD_SOLUTION);
if isempty(fwds)
    fclose(fid);
    error(me,'No forward solutions in %s',fname);
end
%
%   Parent MRI data
%
parent_mri = fiff_dir_tree_find(tree,FIFF.FIFFB_MNE_PARENT_MRI_FILE);
if isempty(parent_mri)
    fclose(fid);
    error(me,'No parent MRI information in %s',fname);
end
%
%   Parent MEG data
%
parent_meg = fiff_dir_tree_find(tree,FIFF.FIFFB_MNE_PARENT_MEAS_FILE);
if isempty(parent_meg)
    fclose(fid);
    error(me,'No parent MEG information in %s',fname);
end
p = 0;
for k = 1:parent_meg.nent
    kind = parent_meg.dir(k).kind;
    pos  = parent_meg.dir(k).pos;
    if kind == FIFF.FIFF_CH_INFO
        p = p+1;
        tag = fiff_read_tag(fid,pos);
        chs(p) = tag.data;
    end
end

try
    src = mne_read_source_spaces(fid,false,tree);
catch
    fclose(fid);
    error(me,'Could not read the source spaces (%s)',mne_omit_first_line(lasterr));
end
for k = 1:length(src)
    src(k).id = mne_find_source_space_hemi(src(k));
end
fwd = [];
%
%   Locate and read the forward solutions
%
megnode = [];
eegnode = [];
for k = 1:length(fwds)
    tag = find_tag(fwds(k),FIFF.FIFF_MNE_INCLUDED_METHODS);
    if isempty(tag)
        fclose(fid);
        error(me,'Methods not listed for one of the forward solutions');
    end
    if tag.data == FIFF.FIFFV_MNE_MEG
        megnode = fwds(k);
    elseif tag.data == FIFF.FIFFV_MNE_EEG
        eegnode = fwds(k);
    end
end

megfwd = read_one(megnode);
if ~isempty(megfwd)
    if megfwd.source_ori == FIFF.FIFFV_MNE_FIXED_ORI
        ori = 'fixed';
    else
        ori = 'free';
    end
    fprintf(1,'\tRead MEG forward solution (%d sources, %d channels, %s orientations)\n',...
        megfwd.nsource,megfwd.nchan,ori);
end
eegfwd = read_one(eegnode);
if ~isempty(eegfwd)
    if eegfwd.source_ori == FIFF.FIFFV_MNE_FIXED_ORI
        ori = 'fixed';
    else
        ori = 'free';
    end
    fprintf(1,'\tRead EEG forward solution (%d sources, %d channels, %s orientations)\n',...
        eegfwd.nsource,eegfwd.nchan,ori);
end
%
%   Merge the MEG and EEG solutions together
%
if ~isempty(megfwd) && ~isempty(eegfwd)
    if size(megfwd.sol.data,2) ~= size(eegfwd.sol.data,2) || ...
            megfwd.source_ori ~= eegfwd.source_ori || ...
            megfwd.nsource ~= eegfwd.nsource || ...
            megfwd.coord_frame ~= eegfwd.coord_frame
        fclose(fid);
        error(me,'The MEG and EEG forward solutions do not match');
    end
    fwd = megfwd;
    fwd.sol.data      = [ fwd.sol.data ; eegfwd.sol.data ];
    fwd.sol.nrow      = fwd.sol.nrow + eegfwd.sol.nrow;
    fwd.sol.row_names = [ fwd.sol.row_names eegfwd.sol.row_names ];
    if  ~isempty(fwd.sol_grad)
        fwd.sol_grad.data      = [ fwd.sol_grad.data ; eegfwd.sol_grad.data ];
        fwd.sol_grad.nrow      = fwd.sol_grad.nrow + eegfwd.sol_grad.nrow;
        fwd.sol_grad.row_names = [ fwd.sol_grad.row_names eegfwd.sol_grad.row_names ];
    end
    fwd.nchan         = fwd.nchan + eegfwd.nchan;
    fprintf(1,'\tMEG and EEG forward solutions combined\n');
elseif ~isempty(megfwd)
    fwd = megfwd;
else
    fwd = eegfwd;
end
clear('megfwd');
clear('eegfwd');
%
%   Get the MRI <-> head coordinate transformation
%
tag = find_tag(parent_mri,FIFF.FIFF_COORD_TRANS);
if isempty(tag)
    fclose(fid);
    error(me,'MRI/head coordinate transformation not found');
else
    mri_head_t = tag.data;
    if mri_head_t.from ~= FIFF.FIFFV_COORD_MRI || mri_head_t.to ~= FIFF.FIFFV_COORD_HEAD
        mri_head_t = fiff_invert_transform(mri_head_t);
        if mri_head_t.from ~= FIFF.FIFFV_COORD_MRI || mri_head_t.to ~= FIFF.FIFFV_COORD_HEAD
            fclose(fid);
            error(me,'MRI/head coordinate transformation not found');
        end
    end
end
fclose(fid);
%
fwd.mri_head_t = mri_head_t;
%
%   Transform the source spaces to the correct coordinate frame
%   if necessary
%
if fwd.coord_frame ~= FIFF.FIFFV_COORD_MRI && ...
        fwd.coord_frame ~= FIFF.FIFFV_COORD_HEAD
    error(me,'Only forward solutions computed in MRI or head coordinates are acceptable');
end
%
nuse = 0;
for k = 1:length(src)
    try
        src(k) = mne_transform_source_space_to(src(k),fwd.coord_frame,mri_head_t);
    catch
        error(me,'Could not transform source space (%s)',mne_omit_first_line(lasterr));
    end
    nuse = nuse + src(k).nuse;
end
if nuse ~= fwd.nsource
    error(me,'Source spaces do not match the forward solution.\n');
end
fprintf(1,'\tSource spaces transformed to the forward solution coordinate frame\n');
fwd.src = src;
%
%   Handle the source locations and orientations
%
if fwd.source_ori == FIFF.FIFFV_MNE_FIXED_ORI || ...
        force_fixed == true
    nuse = 0;
    fwd.source_rr = zeros(fwd.nsource,3);
    fwd.source_nn = zeros(fwd.nsource,3);
    for k = 1:length(src)
        fwd.source_rr(nuse+1:nuse+src(k).nuse,:) = src(k).rr(src(k).vertno,:);
        fwd.source_nn(nuse+1:nuse+src(k).nuse,:) = src(k).nn(src(k).vertno,:);
        nuse = nuse + src(k).nuse;
    end
    %
    %   Modify the forward solution for fixed source orientations
    %
    if fwd.source_ori ~= FIFF.FIFFV_MNE_FIXED_ORI
        fprintf(1,'\tChanging to fixed-orientation forward solution...');
        fix_rot        = mne_block_diag(fwd.source_nn',1);
        fwd.sol.data   = fwd.sol.data*fix_rot;
        fwd.sol.ncol   = fwd.nsource;
        fwd.source_ori = FIFF.FIFFV_MNE_FIXED_ORI;

        if  ~isempty(fwd.sol_grad)
            fwd.sol_grad.data   = fwd.sol_grad.data*kron(fix_rot,eye(3));
            fwd.sol_grad.ncol   = 3*fwd.nsource;
        end
        fprintf(1,'[done]\n');
    end
else
    if surf_ori
        %
        %   Rotate the local source coordinate systems
        %
        fprintf(1,'\tConverting to surface-based source orientations...');
        nuse = 0;
        pp   = 1;
        fwd.source_rr = zeros(fwd.nsource,3);
        for k = 1:length(src)
            fwd.source_rr(nuse+1:nuse+src(k).nuse,:) = src(k).rr(src(k).vertno,:);
            for p = 1:src(k).nuse
                %
                %  Project out the surface normal and compute SVD
                %
                nn = src(k).nn(src(k).vertno(p),:)';
                [ U, S, V ]  = svd(eye(3,3) - nn*nn');
                %
                %  Make sure that ez is in the direction of nn
                %
                if nn'*U(:,3) < 0
                    U = -U;
                end
                fwd.source_nn(pp:pp+2,:) = U';
                pp = pp + 3;
            end
            nuse = nuse + src(k).nuse;
        end
        surf_rot = mne_block_diag(fwd.source_nn',3);
        fwd.sol.data = fwd.sol.data*surf_rot;
        if  ~isempty(fwd.sol_grad)
            fwd.sol_grad.data = fwd.sol_grad.data*kron(surf_rot,eye(3));
        end
        fprintf(1,'[done]\n');
    else
        fprintf(1,'\tCartesian source orientations...');
        nuse = 0;
        fwd.source_rr = zeros(fwd.nsource,3);
        for k = 1:length(src)
            fwd.source_rr(nuse+1:nuse+src(k).nuse,:) = src(k).rr(src(k).vertno,:);
            nuse = nuse + src(k).nuse;
        end
        fwd.source_nn = kron(ones(fwd.nsource,1),eye(3,3));
        fprintf(1,'[done]\n');
    end
end
%
% Add channel information
%
fwd.chs = chs;
%
%   Do the channel selection
%
if ~isempty(include) || ~isempty(exclude)
    %
    %   First do the channels to be included
    %
    if isempty(include)
        pick = ones(1,fwd.nchan);
    else
        pick = zeros(1,fwd.nchan);
        for k = 1:size(include,2)
            c = strmatch(include{k},fwd.sol.row_names,'exact');
            for p = 1:length(c)
                pick(c(p)) = 1;
            end
        end
    end
    %
    %   Then exclude what needs to be excluded
    %
    if ~isempty(exclude)
        for k = 1:length(exclude)
            c = strmatch(exclude{k},fwd.sol.row_names,'exact');
            for p = 1:length(c)
                pick(c(p)) = 0;
            end
        end
    end
    %
    %   Do we have something?
    %
    nuse = sum(pick);
    if nuse == 0
        error(me,'Nothing remains after picking');
    end
    %
    %   Create the selector
    %
    sel = zeros(1,nuse);
    p = 0;
    for k = 1:fwd.nchan
        if pick(k) > 0
            p = p + 1;
            sel(p) = k;
        end
    end
    fprintf(1,'\t%d out of %d channels remain after picking\n',nuse,fwd.nchan);
    %
    %   Pick the correct rows of the forward operator
    %
    fwd.nchan         = nuse;
    fwd.sol.data      = fwd.sol.data(sel,:);
    fwd.sol.nrow      = nuse;
    fwd.sol.row_names = fwd.sol.row_names(sel);

    if ~isempty(fwd.sol_grad)
        fwd.sol_grad.data      = fwd.sol_grad.data(sel,:);
        fwd.sol_grad.nrow      = nuse;
        fwd.sol_grad.row_names = fwd.sol_grad.row_names(sel);
    end

    fwd.chs = fwd.chs(sel);
end

return;

    function [tag] = find_tag(node,findkind)

        for p = 1:node.nent
            if node.dir(p).kind == findkind
                tag = fiff_read_tag(fid,node.dir(p).pos);
                return;
            end
        end
        tag = [];
        return;
    end

    function [one] = read_one(node)
        %
        %   Read all interesting stuff for one forward solution
        %
        if isempty(node)
            one = [];
            return;
        end
        tag = find_tag(node,FIFF.FIFF_MNE_SOURCE_ORIENTATION);
        if isempty(tag)
            fclose(fid);
            error(me,'Source orientation tag not found');
        end
        one.source_ori = tag.data;
        tag = find_tag(node,FIFF.FIFF_MNE_COORD_FRAME);
        if isempty(tag)
            fclose(fid);
            error(me,'Coordinate frame tag not found');
        end
        one.coord_frame = tag.data;
        tag = find_tag(node,FIFF.FIFF_MNE_SOURCE_SPACE_NPOINTS);
        if isempty(tag)
            fclose(fid);
            error(me,'Number of sources not found');
        end
        one.nsource = tag.data;
        tag = find_tag(node,FIFF.FIFF_NCHAN);
        if isempty(tag)
            fclose(fid);
            error(me,'Number of channels not found');
        end
        one.nchan = tag.data;
        try
            one.sol = mne_transpose_named_matrix(fiff_read_named_matrix(fid,node,FIFF.FIFF_MNE_FORWARD_SOLUTION));
        catch
            fclose(fid);
            error(me,'Forward solution data not found (%s)',mne_omit_first_line(lasterr));
        end
        try
            one.sol_grad = mne_transpose_named_matrix(fiff_read_named_matrix(fid,node,FIFF.FIFF_MNE_FORWARD_SOLUTION_GRAD));
        catch
            one.sol_grad = [];
        end
        if size(one.sol.data,1) ~= one.nchan || ...
                (size(one.sol.data,2) ~= one.nsource && size(one.sol.data,2) ~= 3*one.nsource)
            fclose(fid);
            error(me,'Forward solution matrix has wrong dimensions');
        end
        if ~isempty(one.sol_grad)
            if size(one.sol_grad.data,1) ~= one.nchan || ...
                    (size(one.sol_grad.data,2) ~= 3*one.nsource && size(one.sol_grad.data,2) ~= 3*3*one.nsource)
                fclose(fid);
                error(me,'Forward solution gradient matrix has wrong dimensions');
            end
        end
        return;
    end

end
