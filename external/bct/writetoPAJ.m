function writetoPAJ(CIJ, fname, arcs)
%WRITETOPAJ         Write to Pajek
%
%   writetoPAJ(CIJ, fname, arcs);
%
%   This function writes a Pajek .net file from a MATLAB matrix
%
%   Inputs:     CIJ,        adjacency matrix
%               fname,      filename minus .net extension
%               arcs,       1 for directed network
%                           0 for an undirected network
%
%   Chris Honey, Indiana University, 2007


N = size(CIJ,1);
fid = fopen(cat(2,fname,'.net'), 'w');

%%%VERTICES
fprintf(fid, '*vertices %6i \r', N);
for i = 1:N
    fprintf(fid, '%6i "%6i" \r', [i i]);
end

%%%ARCS/EDGES
if arcs
    fprintf(fid, '*arcs \r');
else
    fprintf(fid, '*edges \r');
end

for i = 1:N
    for j = 1:N
        if CIJ(i,j) ~= 0
            fprintf(fid, '%6i %6i %6f \r', [i j CIJ(i,j)]);
        end
    end
end

fclose(fid);

