function  [R,D] = reachdist(CIJ)
%REACHDIST      Reachability and distance matrices
%
%   [R,D] = reachdist(CIJ);
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
%   Note: faster but more memory intensive than "breadthdist.m".
%
%   Algorithm: algebraic path count.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

% initialize
R = CIJ;
D = CIJ;
powr = 2;
N = size(CIJ,1);
CIJpwr = CIJ;

% Check for vertices that have no incoming or outgoing connections.
% These are "ignored" by 'reachdist'.
id = sum(CIJ,1);       % indegree = column sum of CIJ
od = sum(CIJ,2)';      % outdegree = row sum of CIJ
id_0 = find(id==0);    % nothing goes in, so column(R) will be 0
od_0 = find(od==0);    % nothing comes out, so row(R) will be 0
% Use these columns and rows to check for reachability:
col = setxor(1:N,id_0);
row = setxor(1:N,od_0);

[R,D,powr] = reachdist2(CIJ,CIJpwr,R,D,N,powr,col,row);

% "invert" CIJdist to get distances
D = powr - D+1;

% Put 'Inf' if no path found
D(D==(N+2)) = Inf;
D(:,id_0) = Inf;
D(od_0,:) = Inf;


%----------------------------------------------------------------------------

function  [R,D,powr] = reachdist2(CIJ,CIJpwr,R,D,N,powr,col,row)

% Olaf Sporns, Indiana University, 2002/2008

CIJpwr = CIJpwr*CIJ;
R = double(R | ((CIJpwr)~=0));
D = D+R;

if ((powr<=N)&&(~isempty(nonzeros(R(row,col)==0)))) 
   powr = powr+1;
   [R,D,powr] = reachdist2(CIJ,CIJpwr,R,D,N,powr,col,row); 
end;
