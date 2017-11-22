function saveabaqus(node,face,elem,fname,heading)
%
% saveabaqus(node,fname)
% saveabaqus(node,face,fname)
% saveabaqus(node,face,elem,fname)
%
% save a tetrahedral and/or surface mesh as an ABAQUS input file
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2012/07/05
%
% input:
%      node: input, surface node list, dimension (nn,3)
%      face: input, surface face element list, dimension (be,3)
%      elem: input, tetrahedral element list, dimension (ne,4)
%      fname: output file name
%      heading: optional, a descriptive string for the mesh
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

fid=fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end

fprintf(fid,'*HEADING\n');
if(exist('heading','var'))
    fprintf(fid,'**%s\n',heading);
end
fprintf(fid,'*PREPRINT,MODEL=NO,HISTORY=NO,ECHO=NO\n');

if(~isempty(node))
  node=[(1:size(node,1))' node(:,1:3)];
  fprintf(fid,'*NODE, NSET=MeshNode\n');
  fprintf(fid,'%d,\t%e,\t%e,\t%e\n',node');
end

if(~isempty(elem))
  if(size(elem,2)==4)
     elem(:,5)=0;
     labels=0;
  elseif(size(elem,2)<4)
     error('elem must have at least 4 columns');
  else
     labels=unique(elem(:,5));
  end
  len=length(labels);
  count=0;
  elsetall='';
  for i=1:len
      idx=find(elem(:,5)==labels(i));
      fprintf(fid,'*ELEMENT, ELSET=MeshTetra%d, TYPE=C3D4\n',labels(i));
      fprintf(fid,'%d,\t%d,\t%d,\t%d,\t%d\n',[count+(1:length(idx))' elem(idx,1:4)]');
      count=count+length(idx);
      if(length(elsetall)==0)
          elsetall=['MeshTetra' num2str(labels(i))];
      else
          elsetall=sprintf('%s,MeshTetra%d',elsetall,labels(i));
      end
  end
  if(length(elsetall))
      fprintf(fid,'*ELSET, ELSET=MeshTetraAll\n%s\n',elsetall);
  end
end

if(~isempty(face))
  if(size(face,2)==3)
     face(:,4)=0;
     labels=0;
  elseif(size(face,2)<3)
     error('face must have at least 3 columns');
  else
     labels=unique(face(:,4));
  end
  len=length(labels);
  if(exist('count','var')==0)
     count=0;
  end
  elsetall='';
  for i=1:len
      idx=find(face(:,4)==labels(i));
      fprintf(fid,'*ELEMENT, ELSET=MeshTri%d, TYPE=S3R\n',labels(i));
      fprintf(fid,'%d,\t%d,\t%d,\t%d\n',[count+(1:length(idx))' face(idx,1:3)]');
      count=count+length(idx);
      if(length(elsetall)==0)
          elsetall=['MeshTetra' num2str(labels(i))];
      else
          elsetall=sprintf('%s,MeshTri%d',elsetall,labels(i));
      end
  end
  if(length(elsetall))
      fprintf(fid,'*ELSET, ELSET=MeshTriAll\n%s\n',elsetall);
  end
end

fclose(fid);
