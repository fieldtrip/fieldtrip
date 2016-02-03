function r = pagerank_centrality(A, d, falff)
%PAGERANK_CENTRALITY       PageRank centrality
%
%   r = pagerank_centrality(A, d, falff)
%
%   The PageRank centrality is a variant of eigenvector centrality. This
%   function computes the PageRank centrality of each vertex in a graph.
%
%   Formally, PageRank is defined as the stationary distribution achieved
%   by instantiating a Markov chain on a graph. The PageRank centrality of
%   a given vertex, then, is proportional to the number of steps (or amount
%   of time) spent at that vertex as a result of such a process. 
%
%   The PageRank index gets modified by the addition of a damping factor,
%   d. In terms of a Markov chain, the damping factor specifies the
%   fraction of the time that a random walker will transition to one of its
%   current state's neighbors. The remaining fraction of the time the
%   walker is restarted at a random vertex. A common value for the damping
%   factor is d = 0.85.
%
%   Inputs:     A,      adjacency matrix
%               d,      damping factor
%           falff,      initial page rank probability (non-negative)
%
%   Outputs:    r,      vectors of page rankings
%
%   Note: The algorithm will work well for smaller matrices (number of
%   nodes around 1000 or less) 
%
%   References:
%
%   [1]. GeneRank: Using search engine technology for the analysis of
%   microarray experiments, by Julie L. Morrison, Rainer Breitling, Desmond
%   J. Higham and David R. Gilbert, BMC Bioinformatics, 6:233, 2005.
% 	[2]. Boldi P, Santini M, Vigna S (2009) PageRank: Functional
% 	dependencies. ACM Trans Inf Syst 27, 1-23.
%
%   Xi-Nian Zuo, Institute of Psychology, Chinese Academy of Sciences, 2011
%   Rick Betzel, Indiana University, 2012

N = size(A,1);
if nargin < 3
    norm_falff = ones(N,1)/N;
else
    falff = abs(falff);
    norm_falff = falff/sum(falff);
end

deg = sum(A);
ind = (deg == 0);
deg(ind) = 1;
D1 = zeros(N);
D1(1:(N+1):end) = 1./deg;
B = eye(N) - d*(A*D1);
b = (1-d)*norm_falff;
r = B\b;
r = r/sum(r);