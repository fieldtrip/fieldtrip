function  [CIJ,K] = makefractalCIJ(mx_lvl,E,sz_cl)
%MAKEFRACTALCIJ     Synthetic hierarchical modular network
%
%   [CIJ,K] = makefractalCIJ(mx_lvl,E,sz_cl);
%
%   This function generates a directed network with a hierarchical modular
%   organization. All modules are fully connected and connection density 
%   decays as 1/(E^n), with n = index of hierarchical level.
%
%   Inputs:     mx_lvl,     number of hierarchical levels, N = 2^mx_lvl
%               E,          connection density fall-off per level
%               sz_cl,      size of clusters (power of 2)
%
%   Outputs:    CIJ,        connection matrix
%               K,          number of connections present in the output CIJ
%
%
% Olaf Sporns, Indiana University, 2005/2007

% make a little template
t = ones(2).*2;

% compute N and cluster size
N = 2^mx_lvl;
sz_cl = sz_cl-1;

% n = [0 0 0:mx_lvl-3];

for lvl=1:mx_lvl-1
   CIJ = ones(2^(lvl+1),2^(lvl+1));
   group1 = 1:size(CIJ,1)/2;
   group2 = size(CIJ,1)/2+1:size(CIJ,1);
   CIJ(group1,group1) = t;
   CIJ(group2,group2) = t;
   CIJ = CIJ+ones(size(CIJ,1),size(CIJ,1));
   t = CIJ;
end;
s = size(CIJ,1);
CIJ = CIJ-ones(s,s)-mx_lvl.*eye(s);

% assign connection probablities
ee = mx_lvl-CIJ-sz_cl;
ee = (ee>0).*ee;
prob = (1./(E.^ee)).*(ones(s,s)-eye(s));
CIJ = (prob>rand(N));

% count connections
K = sum(sum(CIJ));

