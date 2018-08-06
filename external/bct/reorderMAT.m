function [MATreordered,MATindices,MATcost] = reorderMAT(MAT,H,cost)
%REORDERMAT         Reorder matrix for visualization
%
%   [MATreordered,MATindices,MATcost] = reorderMAT(MAT,H,cost);
%
%   This function reorders the connectivity matrix in order to place more
%   edges closer to the diagonal. This often helps in displaying community
%   structure, clusters, etc.
%
%   Inputs:     MAT,            connection matrix
%               H,              number of reordering attempts
%               cost,           'line' or 'circ', for shape of lattice
%                               (linear or ring lattice)
%
%               MATreordered    reordered connection matrix
%               MATindices      reordered indices
%               MATcost         cost of reordered matrix
%   
%
%   Olaf Sporns, Indiana University


N = length(MAT);
diagMAT = diag(diag(MAT));
MAT = MAT-diagMAT;

% generate cost function
if strcmp(cost,'line');
    profil = fliplr(normpdf(1:N,0,N/2));
end;
if strcmp(cost,'circ');
    profil = fliplr(normpdf(1:N,N/2,N/4));
end;
COST = toeplitz(profil,profil);

% initialize lowCOST
lowMATcost = sum(sum(COST.*MAT));

% keep track of starting configuration
MATstart = MAT;
starta = 1:N;
   
% reorder
for h=1:H
    a = 1:N;
    % choose two positions at random and flip them
    r = randperm(N);
    a(r(1)) = r(2);
    a(r(2)) = r(1);
    MATcostnew = sum(sum(MAT(a,a).*COST));
    if (MATcostnew < lowMATcost)
        MAT = MAT(a,a);
        r2 = starta(r(2));
        r1 = starta(r(1));
        starta(r(1)) = r2;
        starta(r(2)) = r1;
        lowMATcost = MATcostnew;
    end;
end;	% h

MATreordered = MATstart(starta,starta) + diagMAT(starta,starta);
MATindices = starta;
MATcost = lowMATcost;

