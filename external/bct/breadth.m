function [distance,branch] = breadth(CIJ,source)
%BREADTH        Auxiliary function for breadthdist.m
%
%   [distance,branch] = breadth(CIJ,source);
%
%   Implementation of breadth-first search.
%
%   Input:      CIJ,        binary (directed/undirected) connection matrix
%               source,     source vertex
%
%   Outputs:    distance,   distance between 'source' and i'th vertex
%                           (0 for source vertex)
%               branch,     vertex that precedes i in the breadth-first search tree
%                           (-1 for source vertex)
%        
%   Notes: Breadth-first search tree does not contain all paths (or all 
%   shortest paths), but allows the determination of at least one path with
%   minimum distance. The entire graph is explored, starting from source 
%   vertex 'source'.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

N = size(CIJ,1);

% colors: white, gray, black
white = 0; 
gray = 1; 
black = 2;

% initialize colors
color = zeros(1,N);
% initialize distances
distance = inf*ones(1,N);
% initialize branches
branch = zeros(1,N);

% start on vertex 'source'
color(source) = gray;
distance(source) = 0;
branch(source) = -1;
Q = source;

% keep going until the entire graph is explored
while ~isempty(Q)
   u = Q(1);
   ns = find(CIJ(u,:));
   for v=ns
% this allows the 'source' distance to itself to be recorded
      if (distance(v)==0)
         distance(v) = distance(u)+1;
      end;
      if (color(v)==white)
         color(v) = gray;
         distance(v) = distance(u)+1;
         branch(v) = u;
         Q = [Q v];                                             %#ok<AGROW>
      end;
   end;
   Q = Q(2:length(Q));
   color(u) = black;
end
