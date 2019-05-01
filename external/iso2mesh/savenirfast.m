function savenirfast(v,f,filestub, nodeseg, nodeprop, proptype)
%
% savenirfast(nirfaststruct,filestub)
%    or
% savenirfast(v,f,filestub, nodeseg, proptype, proptype)
%
% save a tetrahedral or surface mesh and associated properties to NIRFAST format
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%      nirfaststruct: a structure storing the NIRFAST mesh data, type 
%         'help readnirfast' to read more; alternatively one can use:
%      v: input, node list, the first 3 columns are the x/y/z positions,
%         the remaining columns are combined with nodeprop as node-based
%         (optical) parameters
%      f: input, tetrahedral or surface element list, dimension (ne,3)
%      filestub: output file stub, output will include multiple files
%          filestub.node: node file
%          filestub.elem: element file to store the surface or tet mesh
%          filestub.param: parameter file
%          filestub.region: node label file
%      nodeseg: optional, an integer label field to group nodes into
%         segmentations, same length as v, number starting from 0; or empty
%      nodeprop: optional, additional nodal parameters, typically defined
%         as mua (1/mm), musp (1/mm) and refractive index (n)l; row number
%         equals to that of v, column number is user-defined
%      proptype: optional, the type of the node-property. by default it is 
%         'stnd' - for standard properties; one can also define multi-row
%         header using a cell-array.
%
% example:
%    [node,face,elem]=meshabox([0 0 0],[10 10 10],0.3,1);
%    savenirfast(node,elem,'test', [], ones(size(node)), 'user');
%    mymesh=readnirfast('test')
%    plotmesh([mymesh.nodes mymesh.bndvtx], mymesh.elements,'x>5')
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin<2)
    error('you must provide at least 2 inputs');
end

if(nargin==2)
    filestub=f;
    node=v.nodes;
    f=v.elements;
    proptype=v.type;
    if(isfield(v,'region'))
        nodeseg=v.region;
    end
    if(isfield(v,'mua'))
        nodeprop=[v.mua v.mus v.ri];
    end
    if(isfield(v,'bndvtx'))
        isboundary=v.bndvtx;
    end
    v=node;
end

if(size(v,2)>3)
    if(nargin>4)
        nodeprop=[v(:,4:end) nodeprop];
    else
        nodeprop=v(:,4:end);
    end
else
    if(nargin<5)
        nodeprop=[];
    end
end

if(nargin<6)
    proptype='stnd';
end

if(nargin<4 || isempty(nodeseg))
    nodeseg=zeros(size(v,1),1);
end

if(nargin<6)
    proptype='stnd';
end

if(size(f,2)>4)
    f(:,5:end)=[];
end

if(size(v,2)<3)
    error('v must contain at least 3 columns, and f must have at least 4 columns');
end

if(~exist('isboundary','var'))
    face=surfedge(f);
    isboundary=ismember(1:size(v,1), face(:));
end

fid=fopen([filestub,'.node'],'wt');
if(fid==-1)
    error('Saving node file failed, check permission or disk space.');
end
fprintf(fid,'%d\t%.16f\t%.16f\t%.16f\n',[isboundary(:) v(:,1:3)]');
fclose(fid);

if(size(f,2)<2 || size(f,2)>4)
    error('element list f must contain 3 or 4 columns');
end

fid=fopen([filestub,'.elem'],'wt');
if(fid==-1)
    error('Saving elem file failed, check permission or disk space.');
end
fprintf(fid,'%6d\t%6d\t%6d\t%6d\t\n',f');
fclose(fid);

if(~isempty(nodeseg))
    if(numel(nodeseg)~=size(v,1))
        error('nodeseg must have the same length as v');
    end
    fid=fopen([filestub,'.region'],'wt');
    if(fid==-1)
        error('Saving regin file failed, check permission or disk space.');
    end
    fprintf(fid,'%d\n',nodeseg(:));
    fclose(fid);
end

if(~isempty(nodeprop))
    if(size(nodeprop,1)~=size(v,1))
        error('nodeprop must have the same row number as v');
    end
    fid=fopen([filestub,'.param'],'wt');
    if(fid==-1)
        error('Saving param file failed, check permission or disk space.');
    end
    if(iscell(proptype))
        proptype=strjoin(proptype,'\n');
    end
    fprintf(fid,[proptype '\n']);
    fprintf(fid,[repmat('%.16f\t', 1, size(nodeprop,2)) '\n'],nodeprop');
    fclose(fid);
end
