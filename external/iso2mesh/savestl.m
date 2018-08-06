function savestl(node,elem,fname,solidname)
%
% savestl(node,elem,fname,solidname)
%
% save a tetrahedral mesh to an STL (Standard Tessellation Language) file
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2012/02/14
%
% input:
%      node: input, surface node list, dimension Nx3
%      elem: input, tetrahedral element list; if size(elem,2)==3, it is a surface
%      fname: output file name
%      solidname: an optional string for the name of the object
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

fid=fopen(fname,'wt');
if(fid==-1)
    error('You do not have permission to save mesh files.');
end

if(nargin<4) 
    solidname='';
end
fprintf(fid,'solid %s\n', solidname);

if(nargin<3)
    error('incomplete input');
end
if(isempty(node) || size(node,2)<3)
    error('invalid node input');
end
if(~isempty(elem))
  if(size(elem,2)>=5)
	elem(:,5:end)=[];
  end
  if(size(elem,2)==4)
	elem=meshreorient(node,elem);
	elem=volface(elem);
  end
  ev=surfplane(node,elem);
  ev=ev(:,1:3)./repmat(sqrt(sum(ev(:,1:3).*ev(:,1:3),2)),1,3);
  len=size(elem,1);
  for i=1:len
      fprintf(fid,['facet normal %e %e %e\n' ...
                   ' outer loop\n' ...
                   '  vertex %e %e %e\n' ...
                   '  vertex %e %e %e\n' ...
                   '  vertex %e %e %e\n' ...
                   ' endloop\n' ...
                   'endfacet\n'], [ev(i,:)',node(elem(i,:),1:3)']);
  end
end

fprintf(fid,'endsolid %s\n', solidname);

fclose(fid);
