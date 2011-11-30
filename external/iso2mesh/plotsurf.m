function hm=plotsurf(node,face,varargin)
%
% hm=plotsurf(node,face,opt)
%
% plot 3D surface meshes
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%      node: node coordinates, dimension (nn,3); if node has a 
%            4th column, it will be used to set the color at each node.
%      face: triangular surface face list; if face has a 4th column,
%            it will be used to separate the surface into 
%            sub-surfaces and display them in different colors.
%      opt:  additional options for the plotting, see plotmesh
%
% output:
%   hm: handle or handles (vector) to the plotted surfaces
%
% example:
%
%   h=plotsurf(node,face);
%   h=plotsurf(node,face,'facecolor','r');
% 
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

rngstate = rand ('state');

if(nargin>=2)
    if(size(face,2)==4)
        randseed=hex2dec('623F9A9E'); % "U+623F U+9A9E"
	if(isoctavemesh) randseed=randseed+3; end
        if(~isempty(getvarfrom('base','ISO2MESH_RANDSEED')))
                randseed=getvarfrom('base','ISO2MESH_RANDSEED');
        end
        rand('state',randseed);

        tag=face(:,4);
		types=unique(tag);
        hold on;
		h=[];
	for i=1:length(types)
	    if(size(node,2)==3)
        	h=[h plotasurf(node,face(find(tag==types(i)),1:3),'facecolor',rand(3,1),varargin{:})];
	    else
                h=[h plotasurf(node,face(find(tag==types(i)),1:3),varargin{:})];
	    end
        end
    else
        h=plotasurf(node,face,varargin{:});
    end
end    
if(~isempty(h)) 
  axis equal;
  if(all(get(gca,'view')==[0 90]))
      view(3);
  end
end
if(~isempty(h) & nargout>=1)
  hm=h;
end

rand ('state',rngstate);

%-------------------------------------------------------------------------
function hh=plotasurf(node,face,varargin)
if(size(node,2)==4)
	if(isoctavemesh)
	    h=trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3),node(:,4),'edgecolor','k',varargin{:});
	else
	    h=trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3),node(:,4),varargin{:});
	end
else
	if(isoctavemesh)
	    h=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'edgecolor','k',varargin{:});
	else
	    h=trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3),varargin{:});
	end
end
if(exist('h','var')) hh=h; end
