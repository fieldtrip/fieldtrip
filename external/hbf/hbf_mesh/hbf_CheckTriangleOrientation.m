function orientation=hbf_CheckTriangleOrientation(p,e,verbose)
%HBF_CHECKTRIANGLEORIENTATION tests, how triangle mesh is oriented
% orientation=HBF_CHECKTRIANGLEORIENTATION(mesh)
% orientation=HBF_CHECKTRIANGLEORIENTATION(points,elements,verbose)
%   mesh:   hbf struct for triangle mesh
%   points: mesh vertices, [N x 3]
%   elements:   mesh triangle description, [M x 3]
%   verbose (optional): 1 for error message, 0 for quiet operation
%
%   orientation = 1 for CounterClockWise, 2 for ClockWise,
%                 0 or -1 for a problem
%
%   In hbf, the triangle orientation needs to be CounterClockWise. See
%   HBF_CORRECTTRIANGLEORIENTATION
%
% v180615 (c) Matti Stenroos

if nargin==1
    e=p.e;
    p=p.p;
end
if nargin<3, verbose=1;end
%take a test triangle that most likely is not in a folded part of the mesh.
midpoints=TriangleMidpoints(p,e);
[~ ,maxind]=max(midpoints(:,1));
pmax=midpoints(maxind,:);
p1=p(e(maxind,1),:);
p2=p(e(maxind,2),:);
p3=p(e(maxind,3),:);
%find a point that is just inside this triangle, regardless of mesh
%orientation
normal=cross(p2-p1,p3-p1);
unormal=normal/norm(normal);
if unormal(1)<0
    unormal=-unormal;
end
sl=max([norm(p2-p1),norm(p3-p1)]);
fp=pmax-unormal*.01*sl;

%compute total solid angle spanned by the mesh in this point that is just
%inside the mesh
omegasum=sum(hbf_SolidAngles(e,p,fp));
%omegasum == -4pi
if abs(omegasum+4*pi)<1e-6 && omegasum<0
    orientation=1;
%omegasum == 4pi    
elseif abs(omegasum-4*pi)<1e-6 && omegasum>0
    orientation=2;
%|omegasum| == 0    
elseif abs(omegasum)<1e-6
    if verbose
        fprintf('\nTriangle testpoint outside mesh.\nThe mesh may have very thin regions or sharp corners, please check manually.\n')
    end
    orientation=0;
else 
    if verbose
        fprintf('\nSolid angle not 4pi or 0 --- mesh open or somehow strange. Please check manually.\n')
    end
    orientation=-1;
end

function midpoints=TriangleMidpoints(nodes,elements)
% Calculates midpoints of the mesh triangles.
p1=nodes(elements(:,1),:);
p2=nodes(elements(:,2),:);
p3=nodes(elements(:,3),:);
midpoints=(p1+p2+p3)/3;
