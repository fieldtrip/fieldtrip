function [mesh,success]=hbf_CorrectTriangleOrientation(mesh)
%HBF_CORRECTTRIANGLEORIENTATION checks orientation of triangles and
%   corrects it to CCW (works only in the simplest case).
%
%[mesh_out,success]=HBF_CORRECTTRIANGLEORIENTATION(mesh_in)
%   mesh_in: the mesh to be checked/corrected; hbf struct for triangle mesh
%
%   mesh_out: the checked/corrected mesh; hbf struct for triangle mesh
%   success:    1 for success, 0 for failure
%
% v160229 (c) Matti Stenroos

orientation=hbf_CheckTriangleOrientation(mesh.p,mesh.e);
%if normals point to wrong direction, just flip them...
if orientation==1
    fprintf('Triangle orientation OK.\n');success=1;
elseif orientation<1
    fprintf('Test failed: triangle orientation arbitrary or mesh problematic. \n');
    success=0;
elseif orientation==2
    fprintf('wrong triangle orientation... ');
    enew=mesh.e(:,[1 3 2]);
    orientation=hbf_CheckTriangleOrientation(mesh.p,enew);
    if orientation==1
        fprintf('flipped triangles, OK.\n');
        mesh.e=enew;
        success=1;
    else
        fprintf('flipping failed.\n');
        success=0;
    end
end

