function D = agreement(ci,buffsz)
%AGREEMENT      agreement matrix from clusters
%
%   D = AGREEMENT(CI) takes as input a set of vertex partitions CI of
%   dimensions [vertex x partition]. Each column in CI contains the
%   assignments of each vertex to a class/community/module. This function
%   aggregates the partitions in CI into a square [vertex x vertex]
%   agreement matrix D, whose elements indicate the number of times any two
%   vertices were assigned to the same class.
%
%   In the case that the number of nodes and partitions in CI is large
%   (greater than ~1000 nodes or greater than ~1000 partitions), the script
%   can be made faster by computing D in pieces. The optional input BUFFSZ
%   determines the size of each piece. Trial and error has found that
%   BUFFSZ ~ 150 works well.
%
%   Inputs,     CI,     set of (possibly) degenerate partitions
%               BUFFSZ, optional second argument to set buffer size
%
%   Outputs:    D,      agreement matrix
%
%   Richard Betzel, Indiana University, 2012
%

%modification history
%09.24.2012 - added loop for big N that makes the function slower but also
% prevents it from maxing out memory.

n = size(ci,2);

if nargin < 2
    buffsz = 1000;
end

if n <= buffsz
    
    ind = dummyvar(ci);
    D = ind*ind';
    
else
    
    a = 1:buffsz:n;
    b = buffsz:buffsz:n;
    
    if length(a) ~= length(b)
        b = [b, n];
    end
    
    x = [a' b'];
    nbuff = size(x,1);
    
    D = zeros(size(ci,1));
    for i = 1:nbuff
       y = ci(:,x(i,1):x(i,2));
       ind = dummyvar(y);
       D = D + ind*ind';
    end
    
end

D = D.*~eye(length(D));
