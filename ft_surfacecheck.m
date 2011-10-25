function ft_surfacecheck(cfg,bnd)
% FT_SURFACECHECK computes some quality indexes for the input surface/s.
% The function returns an error if the desired option is not true and if the 'feedback' option is switched 
% on displays the result as a comment string in the workspace. If the bnd argument contains more than 1 
% surface the function automatically checks the desired property on all surfaces.
% 
% Use as
%   cfg = [];
%   cfg.isclosed = 'yes';
%   cfg.feedback = 'yes';
%   ft_surfacecheck(cfg,bnd)
% 
% Optional input arguments 
%   feedback          = 'yes', 'no' (default)
%   isclosed          = 'yes', 'no' (default)
%   isstarshaped      = 'yes', 'no' (default)
%   isoutwardoriented = 'yes', 'no' (default)
%   isnested          = 'yes', 'no' (default)
% 
% See also FT_SURFACEEXTRACT, FT_SURFACEREFINE

% Copyright (C) 2011, Cristiano Micheli, Robert Oostenveld
% 
% $Id: $

ft_defaults

feedback             = ft_getopt(cfg, 'feedback', 'no');
chkisclosed          = istrue(ft_getopt(cfg, 'isclosed'));
chkisstarshaped      = istrue(ft_getopt(cfg, 'isstarshaped','no'));
chkisoutwardoriented = istrue(ft_getopt(cfg, 'isoutwardoriented','no'));
chkisnested          = istrue(ft_getopt(cfg, 'isnested','no'));

% loop over the boundaries
numboundaries = numel(bnd);
msg = {};

if numboundaries>1
  [nesting,order] = nestmatrix(bnd);
end

for i = 1:numboundaries
  
  % print index for each input boundary
  cnt = 1;
  
  if ~isfield(bnd(i),'tri') && ~isfield(bnd(i),'pnt')
    error('the surface %d is not a valid structure\n',i)
  else
    pnt = bnd(i).pnt;
    tri = bnd(i).tri;
  end
  
  % let's limit the check to triangulated surfaces 
  if ~istriangulated(pnt,tri)
    error('the surface %d is not a valid triangulation\n',i);
  end
  
  if chkisclosed
    if isclosed(pnt,tri) 
      msg{cnt} = sprintf('the surface %d is closed\n',i);
      cnt = cnt + 1;
    else
      error(sprintf('the surface %d is not closed\n',i));
    end
  end
  
  if chkisstarshaped
    if isstarshaped(pnt,tri) 
      msg{cnt} = sprintf('the surface %d is star shaped\n',i);
      cnt = cnt + 1;
    else
      error(sprintf('the surface %d is not star shaped\n',i));
    end
  end
  
  if chkisoutwardoriented
    if isoutwardoriented(pnt,tri) 
      msg{cnt} = sprintf('the surface %d is outwards oriented\n',i);
      cnt = cnt + 1;
    else
      error(sprintf('the surface %d is not outwards oriented\n',i));
    end
  end
  
  if chkisnested
    % by convention the first boundary is the innermost
    if sum(nesting(i,:))>0 && i<numboundaries
      msg{cnt} = sprintf('the surface %d is nested\n',i);
      cnt = cnt + 1;    
    elseif sum(nesting(i,:))==0 && i==numboundaries
      %do nothing
    else
      error(sprintf('the surface %d is not nested\n',i));
    end
  end
  
end

if ~isequal(feedback, 'yes') && ~isempty(msg)
  for j=1:cnt-1
    fprintf(msg{j});
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

function [nesting,order] = nestmatrix(bnd)
% nestmatrix checks whether a boundary is nested with all other present
% boundaries

% determine the nesting of the compartments
order = [];
numboundaries = numel(bnd);
nesting = zeros(numboundaries);

for i=1:numboundaries
  for j=1:numboundaries
    if i~=j
      % determine for a single vertex on each surface if it is inside or outside the other surfaces
      curpos = bnd(i).pnt(1,:); % any point on the boundary is ok
      curpnt = bnd(j).pnt;
      curtri = bnd(j).tri;
      nesting(i,j) = bounding_mesh(curpos, curpnt, curtri);
    end
  end
end

if sum(nesting(:))~=(numboundaries*(numboundaries-1)/2)
  error('the compartment nesting cannot be determined');
end

% for a three compartment model, the nesting matrix should look like
%    0 1 1     the first is nested inside the 2nd and 3rd, i.e. the inner skull
%    0 0 1     the second is nested inside the 3rd, i.e. the outer skull
%    0 0 0     the third is the most outside, i.e. the skin
[~, order] = sort(-sum(nesting,2));

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
