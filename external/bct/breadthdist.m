function  [R,D] = breadthdist(CIJ)
%BREADTHDIST      Reachability and distance matrices
%
%   [R,D] = breadthdist(CIJ);
%
%   The binary reachability matrix describes reachability between all pairs
%   of nodes. An entry (u,v)=1 means that there exists a path from node u
%   to node v; alternatively (u,v)=0.
%
%   The distance matrix contains lengths of shortest paths between all
%   pairs of nodes. An entry (u,v) represents the length of shortest path 
%   from node u to  node v. The average shortest path length is the 
%   characteristic path length of the network.
%
%   Input:      CIJ,     binary (directed/undirected) connection matrix
%
%   Outputs:    R,       reachability matrix
%               D,       distance matrix
%
%   Note: slower but less memory intensive than "reachdist.m".
%
%   Algorithm: Breadth-first search.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

N = size(CIJ,1);

D = zeros(N);
for i=1:N
   D(i,:) = breadth(CIJ,i);
end;

% replace zeros with 'Inf's
D(D==0) = Inf;

% construct R
R = double(D~=Inf);

