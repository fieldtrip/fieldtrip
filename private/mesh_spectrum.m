% Mesh spectrum

function [L,H,d] = ct_mesh_spectrum(S,n,varargin)

%[L,H,d] = ct_mesh_spectrum(S,n,mode)
% Compute the mesh laplace matrix and its spectrum
% input,
% S: mesh file, it has to have a pnt and a tri field
% n: number of mesh harmonic functions
% mode: 'full' for the full graph, 'half' if you want to do the first and
% the second half independently (this is useful if your graph is composed
% by two connected components)
% output,
% L: mesh laplacian matrix
% H: matrix containing a mesh harmonic functions per column
% d: spectrum of the negative Laplacian matrix, its units are 1/space^2
% (spatial frequencies are obtained as sqrt(d))


if nargin==2||varargin{1}==1
    pnt{1} = S.pos;
    tri{1} = S.tri;
elseif varargin{1}==2
    pnt{1} = S.pos(1:end/2,:);
    tri{1} = S.tri(1:end/2,:);
    pnt{2} = S.pos(end/2+1:end,:);
    tri{2} = S.tri(end/2+1:end,:) - size(pnt{1},1);
end
   

for j = 1:length(pnt)
    if length(pnt)==2&&j == 1
        disp('Computing the spectrum of the the first hemisphere')
    elseif length(pnt)==2&&j == 2
        disp('Computing the spectrum of the the second hemisphere')
    end
    [L{j},~] = mesh_laplacian(pnt{j},tri{j});
    L{j} = (L{j} + L{j}')/2;
    disp('Computing the spectrum of the negative Laplacian matrix')
    [H{j},D] = eigs(L{j},n,'sm');
    d{j} = diag(D);
    disp('Diagonalization completed')
end
