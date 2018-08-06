function  [CIJ] = makeevenCIJ(N,K,sz_cl)
%MAKEEVENCIJ        Synthetic modular small-world network
%
%   CIJ = makeevenCIJ(N,K,sz_cl);
%
%   This function generates a random, directed network with a specified 
%   number of fully connected modules linked together by evenly distributed
%   remaining random connections.
%
%   Inputs:     N,      number of vertices (must be power of 2)
%               K,      number of edges
%               sz_cl,  size of clusters (power of 2)
%
%   Outputs:    CIJ,    connection matrix
%
%   Notes:  N must be a power of 2.
%           A warning is generated if all modules contain more edges than K.
%           Cluster size is 2^sz_cl;
%
%
%   Olaf Sporns, Indiana University, 2005/2007

% compute number of hierarchical levels and adjust cluster size
mx_lvl = floor(log2(N));
sz_cl = sz_cl-1;

% make a stupid little template
t = ones(2).*2;

% check N against number of levels
Nlvl = 2^mx_lvl;
if (Nlvl~=N) 
    disp('Warning: N must be a power of 2'); 
end;
N = Nlvl;

% create hierarchical template
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

% assign connection probabilities
%CIJp = mx_lvl-CIJ-sz_cl;
%CIJp = (CIJp>0).*CIJp;
CIJp = (CIJ>=(mx_lvl-sz_cl));

% determine number of remaining (non-cluster) connections and their
% possible positions
%CIJc = (CIJp==0);
CIJc = (CIJp==1);
remK = K-nnz(CIJc);
if (remK<0) 
    disp('Warning: K is too small, output matrix contains clusters only');
end;
[a,b] = find(~(CIJc+eye(N)));

% assign 'remK' randomly distributed connections
rp = randperm(length(a));
a = a(rp(1:remK));
b = b(rp(1:remK));
for i=1:remK
   CIJc(a(i),b(i)) = 1;
end;

% prepare for output
CIJ = CIJc;

