function [surf] = mne_read_bem_surfaces(source,add_geom,tree)
%
% [surf] = mne_read_bem_surfaces(source,add_geom,tree)
%
% Reads source spaces from a fif file
%
% source      - The name of the file or an open file id
% add_geom    - Add geometry information to the surfaces
% tree        - Required if source is an open file id, search for the
%               BEM surfaces here
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%   Revision 1.3  2006/09/27 13:07:25  msh
%   Added vertex normals to BEM surfaces.
%   Fixed some data type casts in mne_read_surfaces
%
%   Revision 1.2  2006/05/03 18:53:06  msh
%   Approaching Matlab 6.5 backward compatibility
%
%   Revision 1.1  2006/04/28 18:01:27  msh
%   Added mne_read_bem_surfaces
%   Fixed errors in mne_read_source_spaces and improved triangle data computation.
%
%

me='MNE:mne_read_bem_surfaces';


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
%   These fiff definitions are not needed elsewhere
%
FIFFB_BEM               = 310;    % BEM data
FIFFB_BEM_SURF          = 311;    % One of the surfaces
%
FIFF_BEM_SURF_ID        = 3101;   % int    surface number
FIFF_BEM_SURF_NAME      = 3102;   % string surface name
FIFF_BEM_SURF_NNODE     = 3103;   % int    # of nodes on a surface
FIFF_BEM_SURF_NTRI      = 3104;   % int    # number of triangles on a surface
FIFF_BEM_SURF_NODES     = 3105;   % float  surface nodes (nnode,3)
FIFF_BEM_SURF_TRIANGLES = 3106;   % int    surface triangles (ntri,3)
FIFF_BEM_SURF_NORMALS   = 3107;   % float  surface node normal unit vectors (nnode,3)
FIFF_BEM_COORD_FRAME    = 3112;   % The coordinate frame of the mode
FIFF_BEM_SIGMA          = 3113;   % Conductivity of a compartment
%
%   Default coordinate frame
%
coord_frame = FIFF.FIFFV_COORD_MRI;
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
%   Find BEM
%
bem = fiff_dir_tree_find(tree,FIFFB_BEM);
if isempty(bem)
    close_file;
    error(me,'BEM data not found');
end
bem = bem(1);
%
%   Locate all surfaces
%
bemsurf = fiff_dir_tree_find(bem,FIFFB_BEM_SURF);
if isempty(bemsurf)
    close_file;
    error(me,'BEM surface data not found');
end
fprintf(1,'\t%d BEM surfaces found\n',length(bemsurf));
%
%   Coordinate frame possibly at the top level
%
tag = find_tag(bem,FIFF_BEM_COORD_FRAME);
if ~isempty(tag)
    coord_frame = tag.data;
end
%
%   Read all surfaces
%
for k = 1:length(bemsurf)
    fprintf(1,'\tReading a surface...');
    this = read_source_space(bemsurf(k),coord_frame);
    fprintf(1,'[done]\n');
    if add_geom
        complete_surface_info();
    end
    surf(k) = this;
end

fprintf(1,'\t%d BEM surfaces read\n',length(surf));

close_file();

return;

    function [res] = read_source_space(this,def_coord_frame)
        %
        %   Read all the interesting stuff
        %
        tag = find_tag(this,FIFF_BEM_SURF_ID);
        if isempty(tag)
            res.id = FIFF.FIFFV_BEM_SURF_ID_UNKNOWN;
        else
            res.id = tag.data;
        end

        tag = find_tag(this,FIFF_BEM_SIGMA);
        if isempty(tag)
            res.sigma = 1.0;
        else
            res.sigma = tag.data;
        end

        tag = find_tag(this,FIFF_BEM_SURF_NNODE);
        if isempty(tag)
            close_file();
            error(me,'Number of vertices not found');
        end
        res.np = tag.data;

        tag = find_tag(this,FIFF_BEM_SURF_NTRI);
        if isempty(tag)
            close_file();
            error(me,'Number of triangles not found');
        else
            res.ntri = tag.data;
        end

        tag = find_tag(this,FIFF.FIFF_MNE_COORD_FRAME);
        if isempty(tag)
            tag = find_tag(this,FIFF_BEM_COORD_FRAME);
            if isempty(tag)
                res.coord_frame = def_coord_frame;
            else
                res.coord_frame = tag.data;
            end
        else
            res.coord_frame = tag.data;
        end
        %
        %   Vertices, normals, and triangles
        %
        tag = find_tag(this,FIFF_BEM_SURF_NODES);
        if isempty(tag)
            close_file();
            error(me,'Vertex data not found');
        end
        res.rr = tag.data;
        if size(res.rr,1) ~= res.np
            close_file();
            error(me,'Vertex information is incorrect');
        end

        tag = find_tag(this,FIFF.FIFF_MNE_SOURCE_SPACE_NORMALS);
        if isempty(tag)
            res.nn = [];
        else
            res.nn = tag.data;
            if size(res.nn,1) ~= res.np
                close_file();
                error(me,'Vertex normal information is incorrect');
            end
        end
        tag = find_tag(this,FIFF_BEM_SURF_TRIANGLES);
        if isempty(tag)
            close_file();
            error(me,'Triangulation not found');
        end
        res.tris = tag.data;
        if size(res.tris,1) ~= res.ntri
            close_file();
            error(me,'Triangulation information is incorrect');
        end
        return;


    end

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

    function complete_surface_info()
        %
        %   Main triangulation
        %
        fprintf(1,'\tCompleting triangulation info...');
        fprintf(1,'triangle normals...');
        this.tri_area = zeros(1,this.ntri);
        r1 = this.rr(this.tris(:,1),:);
        r2 = this.rr(this.tris(:,2),:);
        r3 = this.rr(this.tris(:,3),:);
        this.tri_cent = (r1+r2+r3)/3.0;
        this.tri_nn   = cross((r2-r1),(r3-r1));
        %
        %   Triangle normals and areas
        %
        for p = 1:this.ntri
            size = sqrt(this.tri_nn(p,:)*this.tri_nn(p,:)');
            this.tri_area(p) = size/2.0;
            if size > 0.0
                this.tri_nn(p,:) = this.tri_nn(p,:)/size;
            end
        end
        %
        %   Accumulate the vertex normals
        %
        fprintf(1,'vertex normals...');
        this.nn = zeros(this.np,3);
        for p = 1:this.ntri
            this.nn(this.tris(p,:),:) = this.nn(this.tris(p,:),:) + kron(ones(3,1),this.tri_nn(p,:));
        end
        %
        %   Compute the lengths of the vertex normals and scale
        %
        fprintf(1,'normalize...');
        for p = 1:this.np
            size = sqrt(this.nn(p,:)*this.nn(p,:)');
            if size > 0
                this.nn(p,:) = this.nn(p,:)/size;
            end
        end
        fprintf(1,'[done]\n');
    end

    function close_file()

        if open_here
            fclose(fid);
        end

    end


end
