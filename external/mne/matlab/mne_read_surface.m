function [verts, faces] = mne_read_surface(fname)
%
% [verts, faces] = mne_read_surface(fname)
%
% Reads a FreeSurfer surface file
%
% fname       - The file to read
% verts       - Vertex coordinates in meters
% faces       - The triangle descriptions
%               NOTE: The quad file faces are split by this routine to
%               create triangular meshes for all files
%

%
%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.4  2009/02/03 10:26:12  msh
%   Added reading of 'new' quadrangle files
%
%   Revision 1.3  2006/05/26 03:39:12  msh
%   fread3 was called instead of mne_fread3 for the magic number.
%
%   Revision 1.2  2006/05/22 11:01:47  msh
%   Deleted superfluous text from the comments.
%
%   Revision 1.1  2006/05/22 10:44:44  msh
%   Added surface and curvature reading routines.
%   Fixed error in mne_read_source_spaces: triangle normals were not normalized
%
me='MNE:mne_read_surface';
%
%   The input file is big endian
%
fid = fopen(fname,'rb','ieee-be');

if (fid < 0)
    error(me,'Cannot open file %s', fname);
end;
%
%   Magic numbers to identify QUAD and TRIANGLE files
%
%   QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%   NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;
%
NEW_QUAD_FILE_MAGIC_NUMBER =  16777213;
TRIANGLE_FILE_MAGIC_NUMBER =  16777214;
QUAD_FILE_MAGIC_NUMBER     =  16777215;

magic = mne_fread3(fid) ;
if(magic == QUAD_FILE_MAGIC_NUMBER || magic == NEW_QUAD_FILE_MAGIC_NUMBER)
    nvert = mne_fread3(fid);
    nquad = mne_fread3(fid);
    if magic == QUAD_FILE_MAGIC_NUMBER
        fprintf(1,'\t%s is a quad file (nvert = %d nquad = %d)\n', ...
            fname,nvert,nquad);
    else
        fprintf(1,'\t%s is a new quad file (nvert = %d nquad = %d)\n', ...
            fname,nvert,nquad);
    end
    if magic == QUAD_FILE_MAGIC_NUMBER
        verts = fread(fid, nvert*3, 'int16') ./ 100 ;
    else
        verts = fread(fid, nvert*3, 'single');
    end
    if (nargout > 1)
        quads = fread3_many(fid,nquad*4);
        quads = reshape(quads,4,nquad)';
    end
    %
    %   Face splitting follows
    %
    faces = zeros(2*nquad,3);
    nface = 0;
    for k = 1:nquad
        quad = quads(k,:);
        if quad(1) % 2 == 0
            nface = nface + 1;
            faces(nface,1) = quad(1);
            faces(nface,2) = quad(2);
            faces(nface,3) = quad(4);

            nface = nface + 1;
            faces(nface,1) = quad(3);
            faces(nface,2) = quad(4);
            faces(nface,3) = quad(2);
        else
            nface = nface + 1;
            faces(nface,1) = quad(1);
            faces(nface,2) = quad(2);
            faces(nface,3) = quad(3);

            nface = nface + 1;
            faces(nface,1) = quad(1);
            faces(nface,2) = quad(3);
            faces(nface,3) = quad(4);
        end
    end
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
    s = fgets(fid);
    fgets(fid);
    nvert = fread(fid, 1, 'int32') ;
    nface = fread(fid, 1, 'int32') ;
    fprintf(1,'\t%s is a triangle file (nvert = %d ntri = %d)\n',fname,nvert,nface);
    fprintf(1,'\t%s',s);
    verts = fread(fid, nvert*3, 'float32') ;
    faces = fread(fid, nface*3, 'int32') ;
    faces = reshape(faces, 3, nface)';
else
    fclose(fid);
    error(me,'Bad magic number (%d) in surface file %s',magic,fname);
end
verts = 0.001*reshape(verts, 3, nvert)';
faces = faces + 1;                   % Our numbering starts from one
fclose(fid) ;
fprintf(1,'\tRead a surface with %d vertices and %d triangles from %s\n',nvert,nface,fname);

return;

    function [res] = fread3_many(fid,count)

        all = fread(fid,3*count,'uchar');

        res = zeros(count,1);

        p = 1;
        for k = 1:count
            res(k) = bitshift(all(p), 16) + bitshift(all(p+1),8) + all(p+2);
            p = p + 3;
        end
    end

end
