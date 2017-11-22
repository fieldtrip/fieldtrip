function vol=elemvolume(node,elem,option)
%
% vol=elemvolume(node,elem,option)
%
% calculate the volume for a list of simplexes
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/21
%
% input:
%    node:  node coordinates
%    elem:  element table of a mesh
%    option: if option='signed', the volume is the raw determinant,
%            else, the results will be the absolute values
%
% output:
%    vol:   volume values for all elements
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(size(elem,2)==size(node,2))
	enum=size(elem,1);
	vol=zeros(enum,1);
	acol=ones(3,1);
	for i=1:enum
		e1=det([node(elem(i,:),2),node(elem(i,:),3),acol]);
		e2=det([node(elem(i,:),3),node(elem(i,:),1),acol]);
		e3=det([node(elem(i,:),1),node(elem(i,:),2),acol]);
		vol(i)=sqrt(e1*e1+e2*e2+e3*e3)/2;
	end
	return;
end
dim=size(elem,2);
enum=size(elem,1);
vol=zeros(enum,1);
for i=1:enum
    detmat=[node(elem(i,:),:)';ones(1,dim)];
    vol(i)=det(detmat);
end
if(nargin==3 && strcmp(option,'signed'))
  vol=vol/prod(1:size(node,2));
else
  vol=abs(vol)/prod(1:size(node,2));
end
