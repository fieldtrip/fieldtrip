function savedxf(node,face,elem,fname)
%
% savedxf(node,face,elem,fname)
%
% save a surface mesh to DXF format
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2010/04/25
%
% input:
%      node: input, surface node list, dimension (nn,3)
%      face: input, surface face element list, dimension (be,3)
%      elem: input, tetrahedral element list, dimension (ne,4)
%      fname: output file name
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

fprintf(fid,'0\nSECTION\n2\nHEADER\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n');

if(~isempty(face))
  fprintf(fid,'0\nPOLYLINE\n66\n1\n8\nI2M\n70\n64\n71\n%d\n72\n%d\n',size(node,1),size(face,1));
end

if(~isempty(node))
  node=node(:,1:3);
  fprintf(fid,'0\nVERTEX\n8\nI2M\n10\n%f\n20\n%f\n30\n%f\n70\n192\n',node');
end

if(~isempty(face))
  face=face(:,1:3);
  fprintf(fid,'0\nVERTEX\n8\nI2M\n62\n254\n10\n0.0\n20\n0.0\n30\n0.0\n70\n128\n71\n%d\n72\n%d\n73\n%d\n',face');
end

fprintf(fid,'0\nSEQEND\n0\nENDSEC\n');

if(~isempty(face))
  fprintf(fid,'0\nSECTION\n2\nENTITIES\n0\nINSERT\n8\n1\n2\nMesh\n10\n0.0\n20\n0.0\n30\n0.0\n41\n1.0\n42\n1.0\n43\n1.0\n50\n0.0\n0\nENDSEC\n');
end

fprintf(fid,'0\nEOF');

fclose(fid);
