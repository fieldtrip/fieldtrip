function savemphtxt(node, face, elem, filename)
%
% savemphtxt(node, face, elem, filename)
%
% save tetrahedron mesh to comsol file (.mphtxt)
%
% author: Donghyeon Kim (danielkim<at> gist.ac.kr)
% date: 2011/09/29
%
% input:
%      node: input, node list, dimension (nn,3)
%      face: input, surface face element list with label, dimension (be,4)
%      elem: input, tetrahedron element list with label, dimension (ne,5)
%      filename: input, output file name
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

n_node = size(node,1);
n_face = size(face,1);
n_elem = size(elem,1);

if(min(face(:,4))==0)
    face(:,4)=face(:,4)+1;
end
if(min(elem(:,5))==0)
    elem(:,5)=elem(:,5)+1;
end

fp = fopen(filename,'w');
fprintf(fp,'# Created by iso2mesh (http://iso2mesh.sf.net)\n');
fprintf(fp,'0 1\n1\n5 mesh1\n1\n3 obj\n\n');
fprintf(fp,'0 0 1\n4 Mesh\n2\n3\n%d\n1\n',n_node);

% Write Node information
for i = 1:n_node
    fprintf(fp,'%.16f %.16f %.16f\n',node(i,1),node(i,2),node(i,3));
end
fprintf(fp, '\n2\n\n3 tri\n');

% Write Tri information
fprintf(fp, '\n3\n');
fprintf(fp, '%d\n\n',n_face);
for i = 1:n_face
    fprintf(fp,'%d %d %d\n',face(i,1),face(i,2),face(i,3));
end
fprintf(fp, '\n1\n0\n');
fprintf(fp, '%d\n',n_face);
for i = 1:n_face
    fprintf(fp,'%d\n',face(i,4));
end
fprintf(fp, '\n%d\n',n_face);
for i = 1:n_face
    fprintf(fp,'0 0\n');
end

% Write Tet information
fprintf(fp, '\n\n3 tet\n4\n\n%d\n', n_elem);
for i = 1:n_elem
    fprintf(fp,'%d %d %d %d\n',elem(i,1),elem(i,2),elem(i,3),elem(i,4));
end    
fprintf(fp, '\n4\n0\n%d\n%d\n', n_elem);
for i = 1:n_elem
    fprintf(fp,'%d\n',elem(i,5));
end 
fprintf(fp,'\n0\n');
fclose(fp);

