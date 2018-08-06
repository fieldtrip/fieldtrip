function savejmesh(node,face,elem,fname,varargin)
%
% savejmesh(node,face,elem,fname,opt)
%
% export a mesh to the JMesh format
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2011/10/06
%
% input:
%      node: input, node list, dimension (nn,3)
%      face: input, optional, surface face element list, dimension (be,3)
%      elem: input, tetrahedral element list, dimension (ne,4)
%      fname: output file name
%      opt: additional parameters in the form of 'parameter',value pairs
%           valid parameters include:
%           'Dimension': 0 - a user defined mesh, 2- a 2D mesh, 3- a 3D mesh
%           'Author': a string to set the author of the mesh
%           'MeshTitle': a string to set the title of the mesh
%           'MeshTag': a value as the tag of the mesh data
%           'MeshGroup': a value to group this mesh with other mesh data
%           'Comment': a string as the additional note for the mesh data
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin==2)
   fname=face;
   face=[];
   elem=[];
end

if(nargin==3)
   fname=elem;
   elem=[];
end

if(length(varargin)==1 && ischar(varargin{1}))
   opt=struct('FileName',varargin{1});
else
   opt=varargin2struct(varargin{:});
end

meshdim=jsonopt('Dimension',0,opt);

fid=fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end

mesh.MeshVersion=0.5;
mesh.CreationTime=datestr(now);
mesh.Comment=['Created by iso2mesh ' iso2meshver '(http://iso2mesh.sf.net)'];

if(meshdim==0) % a user-defined mesh
    mesh.MeshNode=node;
    if(~isempty(face))
        mesh.MeshSurf=face;
    end
    if(~isempty(elem))
        mesh.MeshElem=elem;
    end
elseif(meshdim==3) % a 3D mesh
    nd=size(node);
    if(nd(2)<3) error('expecting 3 or more columns in node'); end
    mesh.MeshPoint3D=node(:,1:3);
    if(nd(2)>3)
        mesh.MeshNodeVal=node(:,4:end);
    end
    if(~isempty(face))
        nd=size(face);
        if(nd(2)<3) error('expecting 3 or more columns in face'); end
        mesh.MeshTri=face(:,1:3);
        if(nd(2)>3)
            mesh.MeshTriVal=face(:,4:end);
        end
    end
    if(~isempty(elem))
        nd=size(elem);
        if(nd(2)<4) error('expecting 4 or more columns in elem'); end
        mesh.MeshTetra=elem(:,1:4);
        if(nd(2)>4)
            mesh.MeshTetraVal=elem(:,5:end);
        end
    end
elseif(meshdim==2) % a 2D mesh
    nd=size(node);
    if(nd(2)<2) error('expecting 2 or more columns in node'); end
    mesh.MeshPoint2D=node(:,1:2);
    if(nd(2)>2)
        mesh.MeshNodeVal=node(:,3:end);
    end
    if(~isempty(face))
        nd=size(face);
        if(nd(2)<3) error('expecting 3 or more columns in face'); end
        mesh.MeshTri=face(:,1:3);
        if(nd(2)>3)
            mesh.MeshTriVal=face(:,4:end);
        end
    end
    if(~isempty(elem))
        warning('elem is redundant in a 2D mesh, skip');
    end
else
    error('the specified Dimension is not supported, please remove to save data to a general format');
end

author=jsonopt('Author','',opt);
if(~isempty(author))
    mesh.Author=author;
end

title=jsonopt('MeshTitle','',opt);
if(~isempty(title))
    mesh.MeshTitle=title;
end

tag=jsonopt('MeshTag',[],opt);
if(~isempty(tag))
    mesh.MeshTag=tag;
end

group=jsonopt('MeshGroup',[],opt);
if(~isempty(group))
    mesh.MeshTag=group;
end

comm=jsonopt('Comment','',opt);
if(~isempty(comm))
    mesh.Comment=comm;
end

fprintf(fid,'%s\n',savejson('',mesh,varargin{:}));

fclose(fid);
