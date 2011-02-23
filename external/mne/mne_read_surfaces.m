function [surfs] = mne_read_surfaces(surfname,read_curv,read_left,read_right,subject,subjects_dir,add_info)
%
% [surfs] = mne_read_surfaces(surfname,read_curv,read_left,read_right,subject,subjects_dir,add_info)
%
% Reads FreeSurfer surface files for both hemispheres
% as well as curvatures if requested.
%
% Adds the derived geometry information to the surfaces
%
% surfname     - The name of the surface to read, e.g., 'pial'
% read_curv    - read the curvatures as well
% read_left    - read the left hemisphere (default true)
% read_right   - read the right hemisphere (default true)
% subject      - The name of the subject (defaults to SUBJECT environment
%                variable)
% subjects_dir - The name of the MRI data directory (defaults to
%                SUBJECTS_DIR environment variable)
% add_info     - Add auxilliary information to the surfaces
%                (vertex normals, triangle centroids, triangle normals, triangle
%                areas) (default true)
%
% surfs          - Output structure containing the surfaces
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.8  2006/09/27 13:07:25  msh
%   Added vertex normals to BEM surfaces.
%   Fixed some data type casts in mne_read_surfaces
%
%   Revision 1.7  2006/09/25 19:48:16  msh
%   Added projection item kinds to fiff_define_constants
%   Changed some fields to int32 in surface structures
%
%   Revision 1.6  2006/09/15 16:03:02  msh
%   Added an extra argument to mne_read_surfaces to avoid generation of auxilliary data
%
%   Revision 1.5  2006/05/30 18:29:21  msh
%   Return an empty variable if neither read_left nor read_right is set
%
%   Revision 1.4  2006/05/30 18:24:56  msh
%   Added options to read only one hemisphere.
%
%   Revision 1.3  2006/05/22 11:01:47  msh
%   Deleted superfluous text from the comments.
%
%   Revision 1.2  2006/05/22 10:55:02  msh
%   Fixed help text.
%
%   Revision 1.1  2006/05/22 10:44:44  msh
%   Added surface and curvature reading routines.
%   Fixed error in mne_read_source_spaces: triangle normals were not normalized
%

me='MNE:mne_read_surfaces';

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if nargin < 7
    add_info = true;
    if nargin < 6
        subjects_dir=getenv('SUBJECTS_DIR');
        if isempty(subjects_dir)
            error(me,'SUBJECTS_DIR not set');
        end
        if nargin < 5
            subject=getenv('SUBJECT');
            if isempty(subject)
                error(me,'SUBJECT not set');
            end
        end
        if nargin < 4
            read_right = true;
        end
        if nargin < 3
            read_left = true;
        end
    end
elseif nargin ~= 7
    error(me,'Incorrect number of arguments');
end
%
%   Read both hemispheres
%
if read_left
    [ lhvert, lhtri ] = mne_read_surface(sprintf('%s/%s/surf/lh.%s',subjects_dir,subject,surfname));
    fprintf(1,'\n');
end
if read_right
    [ rhvert, rhtri ] = mne_read_surface(sprintf('%s/%s/surf/rh.%s',subjects_dir,subject,surfname));
    fprintf(1,'\n');
end
%
if read_curv
    if read_left
        lhcurv = mne_read_curvature(sprintf('%s/%s/surf/lh.curv', ...
            subjects_dir,subject));
    end
    if read_right
        rhcurv = mne_read_curvature(sprintf('%s/%s/surf/rh.curv', ...
            subjects_dir,subject));
    end
    fprintf(1,'\n');
else
    lhcurv = [];
    rhcurv = [];
end
%
%   Compose the source space structure
%
nsurf = 0;
if read_left
    this.id           = int32(FIFF.FIFFV_MNE_SURF_LEFT_HEMI);
    this.np           = int32(size(lhvert,1));
    this.ntri         = int32(size(lhtri,1));
    this.coord_frame  = int32(FIFF.FIFFV_COORD_MRI);
    this.rr           = lhvert;
    this.nn           = [];
    this.tris         = cast(lhtri,'int32');
    this.nuse         = this.np;
    this.inuse        = ones(1,this.np);
    this.vertno       = [ 1 : this.np ];
    this.curv         = lhcurv;
    %
    %   Add the derived geometry data
    %
    hemi='left';
    if add_info == true
        complete_surface_info;
    end
    nsurf = nsurf + 1;
    surfs(nsurf) = this;
    clear('this');
end
if read_right
    this.id           = int32(FIFF.FIFFV_MNE_SURF_RIGHT_HEMI);
    this.np           = int32(size(rhvert,1));
    this.ntri         = int32(size(rhtri,1));
    this.coord_frame  = int32(FIFF.FIFFV_COORD_MRI);
    this.rr           = rhvert;
    this.nn           = [];
    this.tris         = cast(rhtri,'int32');
    this.nuse         = int32(this.np);
    this.inuse        = int32(ones(1,this.np));
    this.vertno       = int32([ 1 : this.np ]);
    this.curv         = rhcurv;
    hemi='right';
    if add_info == true
        complete_surface_info;
    end
    nsurf = nsurf + 1;
    surfs(nsurf) = this;
    clear('this');
end

if nsurf == 0
    surfs = [];
end

return;

    function complete_surface_info
        %
        %   Main triangulation
        %
        fprintf(1,'\tCompleting triangulation info for the %s hemisphere...',hemi);
        this.tri_area = zeros(1,this.ntri);
        fprintf(1,'triangle normals...');
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
end
