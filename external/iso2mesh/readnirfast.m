function nirfastmesh=readnirfast(filestub)
%
% nirfastmesh=readnirfast(v,f,filestub)
%
% load a group of mesh files saved in the NIRFAST format
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%      filestub: output file stub, output will include multiple files
%          filestub.node: node file
%          filestub.elem: element file to store the surface or tet mesh
%          filestub.param: parameter file
%          filestub.region: node label file
%          filestub.excoef: extinction coeff list
%
% output:
%      nirfastmesh.nodes: node list, 3 columns
%      nirfastmesh.elements: element list, 3 or 4 columns integers
%      nirfastmesh.bndvtx: boundary flag for each node, 1: on the boundary
%      nirfastmesh.region: node segmentation labels
%      nirfastmesh.dimension: dimension of the mesh
%      nirfastmesh.excoef: extinction coeff list
%      nirfastmesh.excoefheader: extinction coeff list field names
%      nirfastmesh.type: the header of the .param file
%      nirfastmesh.prop: optical property list (non-standard, need further processing)
%
%   format definition see http://www.dartmouth.edu/~nir/nirfast/tutorials/NIRFAST-Intro.pdf
%
% example:
%    [node,face,elem]=meshabox([0 0 0],[10 10 10],0.3,1);
%    savenirfast(node,elem,'test', [], ones(size(node)), 'user');
%    mymesh=readnirfast('test')
%    plotmesh([mymesh.nodes mymesh.bndvtx], mymesh.elements,'x>5')
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fname=[filestub,'.node'];
if(~exist(fullfile(pwd,fname),'file'))
    error([fname ' could not be found']);
end
nirfastmesh.nodes=load(fname);

nirfastmesh.bndvtx=nirfastmesh.nodes(:,1);
nirfastmesh.nodes(:,1)=[];

fname=[filestub,'.elem'];
if(~exist(fullfile(pwd,fname),'file'))
    error([fname ' could not be found']);
end
nirfastmesh.elements=load(fname);
nirfastmesh.dimension=size(nirfastmesh.elements,2)-1;

fname=[filestub,'.region'];
if(exist(fullfile(pwd,fname),'file'))
    nirfastmesh.region=load(fname);
end

fname=[filestub,'.excoef'];
fid=fopen(fname,'rt');
if(fid>=0)
    linenum=0;
    textheader={};
    while(~feof(fid))
        oneline=fgetl(fid);
        linenum=linenum+1;
        [data, count]=sscanf(oneline,'%f');
        if(count>1)
            params=fscanf(fid,repmat('%f ',1,count),inf);
            params=reshape(params,length(params)/count, count);
            params(2:end+1,:)=params;
            params(1,:)=data(:)';
            nirfastmesh.excoef=params;
            nirfastmesh.excoefheader=textheader;
            break;
        else
            textheader{end+1}=oneline;
        end
    end
    fclose(fid);
end

fname=[filestub,'.param'];
fid=fopen(fname,'rt');
if(fid>=0)
    linenum=0;
    params=[];
    while(~feof(fid))
        oneline=fgetl(fid);
        if(linenum==0)
            nirfastmesh.type=oneline;
        end
        linenum=linenum+1;
        [data, count]=sscanf(oneline,'%f');
        if(count>1)
            params=fscanf(fid,repmat('%f ',1,count),inf);
            params=reshape(params,length(params)/count, count);
            params(2:end+1,:)=params;
            params(1,:)=data(:)';
            nirfastmesh.prop=params;
            break;
        end
    end
    fclose(fid);
end
