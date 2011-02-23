function [src] = mne_read_source_spaces(source,add_geom,tree)
%
% [src] = mne_read_source_spaces(source,add_geom,tree)
%
% Reads source spaces from a fif file
%
% source      - The name of the file or an open file id
% add_geom    - Add geometry information to the source spaces
% tree        - Required if source is an open file id, search for the
%               source spaces here
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.13  2007/03/11 13:32:22  msh
%   Nearest numbers had to be incremented by one to comply with the 1-based Matlab indexing
%
%   Revision 1.12  2006/09/25 19:48:16  msh
%   Added projection item kinds to fiff_define_constants
%   Changed some fields to int32 in surface structures
%
%   Revision 1.11  2006/05/22 10:44:44  msh
%   Added surface and curvature reading routines.
%   Fixed error in mne_read_source_spaces: triangle normals were not normalized
%
%   Revision 1.10  2006/05/05 03:50:40  msh
%   Added routines to compute L2-norm inverse solutions.
%   Added mne_write_inverse_sol_stc to write them in stc files
%   Several bug fixes in other files
%
%   Revision 1.9  2006/05/03 19:05:06  msh
%   Eliminated more Matlab 6.5 incompatibility
%
%   Revision 1.8  2006/05/03 18:53:06  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.7  2006/05/03 15:50:39  msh
%   Fixed some Matlab 6.5 incompatibility.
%
%   Revision 1.6  2006/04/28 18:01:27  msh
%   Added mne_read_bem_surfaces
%   Fixed errors in mne_read_source_spaces and improved triangle data computation.
%
%   Revision 1.5  2006/04/23 15:29:41  msh
%   Added MGH to the copyright
%
%   Revision 1.4  2006/04/20 21:49:38  msh
%   Added mne_read_inverse_operator
%   Changed some of the routines accordingly for more flexibility.
%
%   Revision 1.3  2006/04/18 23:21:22  msh
%   Added mne_transform_source_space_to and mne_transpose_named_matrix
%
%   Revision 1.2  2006/04/18 20:44:46  msh
%   Added reading of forward solution.
%   Use length instead of size when appropriate
%
%   Revision 1.1  2006/04/17 21:05:08  msh
%   Added source space reading routine.
%
%

me='MNE:mne_read_source_spaces';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin == 1
    add_geom = false;
    fid      = -1;
    fname    = source;
elseif nargin == 2
    fid      = -1;
    fname    = source;
elseif nargin == 3
    a = whos('source');
    if strcmp(a.class,'char')
        fname = source;
        fid   = -1;
    else
        fid   = source;
    end
else
    error(me,'Incorrect number of arguments');
end
%
%   Open the file, create directory
%
if fid < 0
    [ fid, tree ] = fiff_open(fname);
    open_here = true;
else
    open_here = false;
end
%
%   Find all source spaces
%
spaces = fiff_dir_tree_find(tree,FIFF.FIFFB_MNE_SOURCE_SPACE);
if isempty(spaces)
    close_file;
    error(me,'No source spaces found');
end

for k = 1:length(spaces)
    fprintf(1,'\tReading a source space...');
    this = read_source_space(spaces(k));
    fprintf(1,'[done]\n');
    if add_geom
        complete_source_space_info;
    end
    src(k) = this;
end

fprintf(1,'\t%d source spaces read\n',length(spaces));

close_file;

return;

    function [res] = read_source_space(this)
        %
        %   Read all the interesting stuff
        %
        FIFF_BEM_SURF_NTRI=3104;
        FIFF_BEM_SURF_TRIANGLES=3106;
        
        tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_ID);
        if isempty(tag)
            res.id = int32(FIFF.FIFFV_MNE_SURF_UNKNOWN);
        else
            res.id = int32(tag.data);
        end
        
        tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NPOINTS);
        if isempty(tag)
            close_file;
            error(me,'Number of vertices not found');
        end
        res.np = tag.data;
        
        tag = find_tag(this,FIFF_BEM_SURF_NTRI);
        if isempty(tag)
            tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NTRI);
            if isempty(tag)
                res.ntri = 0;
            else
                res.ntri = tag.data;
            end
        else
            res.ntri = tag.data;
        end
        
        tag = find_tag(this,FIFF.FIFF_MNE_COORD_FRAME);
        if isempty(tag)
            close_file;
            error(me,'Coordinate frame information not found');
        end
        res.coord_frame = tag.data;
        %
        %   Vertices, normals, and triangles
        %
        tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_POINTS);
        if isempty(tag)
            close_file;
            error(me,'Vertex data not found');
        end
        res.rr = tag.data;
        if size(res.rr,1) ~= res.np
            close_file;
            error(me,'Vertex information is incorrect');
        end
        
        tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NORMALS);
        if isempty(tag)
            close_file;
            error(me,'Vertex normals not found');
        end
        res.nn = tag.data;
        if size(res.nn,1) ~= res.np
            close_file;
            error(me,'Vertex normal information is incorrect');
        end
        
        if res.ntri > 0
            tag = find_tag(this,FIFF_BEM_SURF_TRIANGLES);
            if isempty(tag)
                tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_TRIANGLES);
                if isempty(tag)
                    close_file;
                    error(me,'Triangulation not found');
                else
                    res.tris = tag.data;
                end
            else
                res.tris = tag.data;
            end
            if size(res.tris,1) ~= res.ntri
                close_file;
                error(me,'Triangulation information is incorrect');
            end
        else
            res.tris = [];
        end
        %
        %   Which vertices are active
        %
        tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NUSE);
        if isempty(tag)
            res.nuse   = 0;
            res.inuse  = int32(zeros(1,res.nuse));
            res.vertno = [];
        else
            res.nuse = tag.data;
            tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_SELECTION);
            if isempty(tag)
                close_file;
                error(me,'Source selection information missing');
            end
            res.inuse  = int32(tag.data');
            res.vertno = int32(zeros(1,res.nuse));
            if length(res.inuse) ~= res.np
                close_file;
                error(me,'Incorrect number of entries in source space selection');
            end
            pp = 0;
            for p = 1:res.np
                if res.inuse(p)
                    pp = pp + 1;
                    res.vertno(pp) = p;
                end
            end
        end
        %
        %   Use triangulation
        %
        tag1 = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NUSE_TRI);
        tag2 = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_USE_TRIANGLES);
        if isempty(tag1) || isempty(tag2)
            res.nuse_tri = int32(0);
            res.use_tris = [];
        else
            res.nuse_tri = tag1.data;
            res.use_tris = tag2.data;
        end
        %
        %   Patch-related information
        %
        tag1 = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NEAREST);
        tag2 = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NEAREST_DIST);
        if isempty(tag1) || isempty(tag2)
            res.nearest = [];
            res.nearest_dist = [];
        else
            res.nearest = tag1.data' + 1;
            res.nearest_dist = tag2.data';
        end
        res.pinfo = mne_patch_info(res.nearest);
        if ~isempty(res.pinfo)
            fprintf(1,'Patch information added...');
        end
        %
        % Distances
        %
        tag1 = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_DIST);
        tag2 = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_DIST_LIMIT);
        if isempty(tag1) || isempty(tag2)
            res.dist = [];
            res.dist_limit = [];
        else
            res.dist       = tag1.data;
            res.dist_limit = tag2.data;
            %
            %  Add the upper triangle
            %
            res.dist = res.dist + res.dist';
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
        
    end

    function complete_source_space_info
        %
        %   Main triangulation
        %
        fprintf(1,'\tCompleting triangulation info...');
        this.tri_area = zeros(1,this.ntri);
        r1 = this.rr(this.tris(:,1),:);
        r2 = this.rr(this.tris(:,2),:);
        r3 = this.rr(this.tris(:,3),:);
        this.tri_cent = (r1+r2+r3)/3.0;
        this.tri_nn   = cross((r2-r1),(r3-r1));
        for p = 1:this.ntri
            size = sqrt(this.tri_nn(p,:)*this.tri_nn(p,:)');
            this.tri_area(p) = size/2.0;
            this.tri_nn(p,:) = this.tri_nn(p,:)/size;
        end
        fprintf(1,'[done]\n');
        %
        %   Selected triangles
        %
        fprintf(1,'\tCompleting selection triangulation info...');
        if this.nuse_tri > 0
            r1 = this.rr(this.use_tris(:,1),:);
            r2 = this.rr(this.use_tris(:,2),:);
            r3 = this.rr(this.use_tris(:,3),:);
            this.use_tri_cent = (r1+r2+r3)/3.0;
            this.use_tri_nn   = cross((r2-r1),(r3-r1));
            for p = 1:this.nuse_tri
                this.use_tri_area(p)   = sqrt(this.use_tri_nn(p,:)*this.use_tri_nn(p,:)')/2.0;
            end
        end
        fprintf(1,'[done]\n');
    end

    function close_file
        
        if open_here
            fclose(fid);
        end
        
    end


end
