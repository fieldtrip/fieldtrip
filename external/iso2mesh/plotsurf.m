function hm=plotsurf(node,face,varargin)
%
% hm=plotsurf(node,face,opt)
%
% plot 3D surface meshes (2d manifold) or polylines (1d manifold)
% 
% author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%
% input: 
%      node: node coordinates, dimension (nn,3); if node has a 
%            4th column, it will be used to set the color at each node.
%      face: triangular surface face list; if face has a 4th column,
%            it will be used to separate the surface into 
%            sub-surfaces and display them in different colors;
%            face can be a cell array, each element of the array represents
%            a polyhedral facet of the mesh, if an element is an array with
%            two array subelements, the first one is the node index, the
%            second one is a scalar as the group id of the facet.
%      opt:  additional options for the plotting, see plotmesh
%
% output:
%   hm: handle or handles (vector) to the plotted surfaces
%
% example:
%
%   h=plotsurf(node,face);
%   h=plotsurf(node,face,'facecolor','r');
%   h=plotsurf(node,edges,'linestyle','-','linewidth',2,'color','r');
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

rngstate = rand ('state');

if(nargin>=2)
  randseed=hex2dec('623F9A9E'); % "U+623F U+9A9E"
  if(isoctavemesh) randseed=randseed+3; end
  if(~isempty(getvarfrom({'caller','base'},'ISO2MESH_RANDSEED')))
        randseed=getvarfrom({'caller','base'},'ISO2MESH_RANDSEED');
  end
  rand('state',randseed);

  if(iscell(face))
    sc=sparse(10,3); % face colormap
    sc(1:10,:)=rand(3,10)';
    len=length(face);
    newsurf=cell(1);
    % reorganizing each labeled surface into a new cell
    for i=1:len
        fc=face{i};
        if(iscell(fc) && length(fc)>=2)
            if(fc{2}+1>10)
                sc(fc{2}+1,:)=rand(1,3);
            end
            if(fc{2}+1>length(newsurf))
                newsurf{fc{2}+1}={};
            end
            newsurf{fc{2}+1}{end+1}=fc{1};
        else % unlabeled facet is tagged by 0
            if(iscell(fc))
                newsurf{1}{end+1}=cell2mat(fc);
            else
                newsurf{1}{end+1}=fc;
            end
        end
    end
    hold on;
	h=[];
    newlen=length(newsurf);

    for i=1:newlen
        if(isempty(newsurf{i})); continue; end
        try 
            subface=cell2mat(newsurf{i}')';
            if(size(subface,1)>1 && ismatrix(subface))
               subface=subface';
            end
            h=[h patch('Vertices',node,'Faces',subface,'facecolor',sc(i,:),varargin{:})];
        catch
            for j=1:length(newsurf{i})
                h=[h patch('Vertices',node,'Faces',newsurf{i}{j},'facecolor',sc(i,:),varargin{:})];
            end
        end
    end
  else
    if(size(face,2)==4)
        tag=face(:,4);
		types=unique(tag);
        hold on;
		h=[];
	    for i=1:length(types)
            if(size(node,2)==3)
                h=[h plotasurf(node,face(tag==types(i),1:3),'facecolor',rand(3,1),varargin{:})];
            else
                h=[h plotasurf(node,face(tag==types(i),1:3),varargin{:})];
            end
        end
    else
        h=plotasurf(node,face,varargin{:});
    end
  end
end    
if(~isempty(h)) 
  axis equal;
  if(all(get(gca,'view')==[0 90]))
      view(3);
  end
end
if(~isempty(h) && nargout>=1)
  hm=h;
end

rand ('state',rngstate);

%-------------------------------------------------------------------------
function hh=plotasurf(node,face,varargin)
isoct=isoctavemesh;
if(size(face,2)<=2)
    h=plotedges(node,face,varargin{:});
else
  if(size(node,2)==4)
	if(isoct && ~exist('trisurf','file'))
	    h=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),node(:,4),'edgecolor','k',varargin{:});
	else
	    h=trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3),node(:,4),varargin{:});
	end
  else
	if(isoct && ~exist('trisurf','file'))
	    h=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'edgecolor','k',varargin{:});
	else
	    h=trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3),varargin{:});
    end
  end
end
if(exist('h','var')) hh=h; end
