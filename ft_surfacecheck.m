function ft_surfacecheck(cfg,bnd)
% FT_SURFACECHECK computes some quality indexes for the input surface
% and displays them in the command line
% 
% Use as
%   ft_surfacecheck(cfg,bnd)
% 
% Optional input arguments 
%   feedback           = yes, no
% 
% See also FT_SURFACEEXTRACT, FT_SURFACEREFINE

% Copyright (C) 2011, Cristiano Micheli, Robert Oostenveld
% 
% $Id: $

ft_defaults

feedback  = ft_getopt(cfg, 'feedback', 'yes');

% let's limit the check to triangulated surfaces 
if ~isfield(bnd,'tri') && ~isfield(bnd,'pnt')
  error('this is not a valid surface structure')
else
  pnt = bnd.pnt;
  tri = bnd.tri;
end

% print some properties
if ~isequal(feedback, 'no')
  if ~istriangulated(pnt,tri)
    fprintf('the surface is not a valid triangulation\n');
    return
  end
  if isclosed(pnt,tri)
    fprintf('the surface is closed\n');
  else
    fprintf('the surface is not closed\n');
  end
  if isstarshaped(pnt,tri)
    fprintf('the surface is star shaped\n');
  else
    fprintf('the surface is not star shaped\n');
  end  
  if isoutwardoriented(pnt,tri)
    fprintf('the surface is outwards oriented\n');
  else
    fprintf('the surface is not outwards oriented\n');
  end
end

function status=istriangulated(pnt,tri)
status = false;
if size(pnt,2)==3 && size(tri,2)==3 && isint(tri)
  status = true;
end

function status=isclosed(pnt,tri)
% IS_CLOSED checks whether a triangulation describes a closed surface
%
% Use as
%    status = isclosed(pnt,tri)
% which will return a single boolean value.
%
% Copyright (C) 2009, Robert Oostenveld

status = false;

npnt = size(pnt,1);
ntri = size(tri,1);

% the number of triangles should be consistent with the number of vertices
if ntri~=2*(npnt-2)
  status = false;
  return
end

% determine a point outside the mesh
ori = max(pnt) + [1 1 1];
% shift that point towards the origin, and the mesh along with it
pnt(:,1) = pnt(:,1) - ori(1);
pnt(:,2) = pnt(:,2) - ori(2);
pnt(:,3) = pnt(:,3) - ori(3);

w = solid_angle(pnt, tri);
w = sum(w);

% the total solid angle should be zero for a closed surface
if abs(w)<100*eps
  status = true;
else
  status = false;
end

function status=isoutwardoriented(pnt,tri)
% ISOUTWARDORIENTED checks whether the triangles in a mesh are all outward oriented
%
% Use as
%   status = isoutward_oriented(pnt,tri)
%
% Copyright (C) 2009, Robert Oostenveld

status = false;
ok = isstarshaped(pnt, tri);

if ~ok
  % shift the geometrical mean of all vertices towards the origin
  ori = mean(pnt,1);
  pnt(:,1) = pnt(:,1) - ori(1);
  pnt(:,2) = pnt(:,2) - ori(2);
  pnt(:,3) = pnt(:,3) - ori(3);
  ok = isstarshaped(pnt, tri);
end

if ~ok
  % shift the center of the bounding box towards the origin
  ori = (min(pnt,1) + max(pnt,1))./2;
  pnt(:,1) = pnt(:,1) - ori(1);
  pnt(:,2) = pnt(:,2) - ori(2);
  pnt(:,3) = pnt(:,3) - ori(3);
  ok = isstarshaped(pnt, tri);
end

if ~ok
  error('cannot determine whether the mesh is inward or outward oriented, the mesh is not starshaped');
else
  w = solid_angle(pnt, tri);
  w = sum(w);

  % the total solid angle of the surface should be either 4*pi or -4*pi
  if w>0 && (abs(w)-4*pi)<1000*eps
    % it is inward oriented
    status = false;
  elseif w<0 && (abs(w)-4*pi)<1000*eps
    % it is outward oriented
    status = true;
  else
    error(sprintf('cannot determine whether the mesh is inward or outward oriented, inconsistend total solid angle (%g)', w));
  end  
end

function status=isstarshaped(pnt,tri,ori,flag)
% IS_STARSHAPED checks whether a triangulation is star-shaped, given the
% center point. If no center point is specified, the center off
% mass of all vertices will be used.
%
% A surface is star-shaped if the complete surface (i.e. all triangles) are visible from the center.
%
% Use as
%   status = isstarshaped(pnt, tri, ori, flag)
% which will return a single boolean value if flag==0 (default) or a boolean vector if flag==1.
% Copyright (C) 2009, Robert Oostenveld
status = false;

if nargin<3 || isempty(ori)
  ori = mean(pnt,1);
end

if nargin<4 || isempty(flag)
  % only return a single boolean value
  flag = 0;
end

% shift the mesh towards the origin
pnt(:,1) = pnt(:,1) - ori(1);
pnt(:,2) = pnt(:,2) - ori(2);
pnt(:,3) = pnt(:,3) - ori(3);

nrm = normals(pnt, tri, 'triangle');

npnt = size(pnt,1);
ntri = size(tri,1);

if flag==0
  status = true;
else
  status = true(ntri,1);
end

for i=1:ntri
% compute the geometrical centre of each triangle
  avg = mean(pnt(tri(i,:),:),1);
  % if this triangle is visible from the center, then the dot-product will be positive
  if dot(avg, nrm(i,:))<0
    if flag==0
      % only return a single boolean value
      status = false;
      return
    else
      % return a boolean value for every triangle
      status(i) = false;
    end
  end
end

function res=isint(v)
% Use as
%   [res, int_part, dec_part] = isint(v)
% Examples:
%   [res, int_part, dec_part] = isint(4.07)
%   [res, int_part, dec_part] = isint(4.00)
%
%     Author : Ziad Saad
%     Date : Wed Aug 11 16:42:55 CDT 1999 

fv = fix(v);
df = v - fv;
if (~df),
	res = 1;
	dec_part = 0;
	int_part = v;
else
	res = 0;
	dec_part = df.*sign(v);
	int_part = fv.*sign(v);
end
