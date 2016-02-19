function [d] = triangle2distance(tri, pos, s)

% TRIANGLE2DISTANCE computes the geodesic distance (across the edges) on a
% mesh, using Dijkstra's algorithm. The Dijkstra code is an efficient
% vectorized version of a function from MIT's graphtool toolbox, operating
% on an adjacency matrix.
%
% Use as
%   d = triangle2distance(tri, pos, s)
%
% Input arguments:
%   tri = Mx3 matrix describing the triangles
%   pos = Nx3 matrix describing the position of the vertices
%   s   = (can be empty), scalar or vector with indices for the points for
%         which the distance (to all other points) will be computed. If
%         empty or not defined, all points will be considered.
%
% Output argument:
%   d   = Nxnumel(s) distance matrix

% Copyright (C) 2015, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id:$

adj = triangle2connectivity(tri, pos);
n   = length(adj);
if nargin==2
  s  = 1:n;
end
ns = length(s);


d    = inf*ones(n,ns); % distance s-all nodes
for k = 1:ns
  d(s(k),k) = 0; % s-s distance
end

for k = 1:ns
  if mod(k,10)==0,
    fprintf('computing distance between node %d and nodes 1:%d\n', s(k), n);
  end
  T = 1:n;    % node set with shortest paths not found
  %while ~isempty(T) %%% a for-loop goes ~10% faster, and we know how often
  %we need to iterate
  for m = 1:n  
    [dmin,ind] = min(d(T,k));
    
    adj_ = adj(T,T(ind));
    d_   = d(T,k);
    
    % logic: shrink the distance if there's an edge (adj_>0) AND if it's
    % shorter to travel through this edge
    criterion = adj_ > 0 & d_ > d(T(ind),k)+adj_;
    d(T(criterion),k) = d(T(ind),k) + adj_(criterion);
    
    T(ind) = [];
  end
end

% the below code is from MIT's graphtool toolbox. the above code that
% computes the distance based on an adjacency matrix is taken from there,
% but the vectorized functionality is ~10 times as fast.


% Implements a simple version of the Dijkstra shortest path algorithm
% Returns the distance from a single vertex to all others, doesn't save the path
% INPUTS: adjacency matrix (adj), start node (s)
% OUTPUTS: shortest path length from start node to all other nodes
% Note: works with a weighted/directed matrix
% GB, Last Updated: December 13, 2004

% Copyright (c) 2011, Massachusetts Institute of Technology.
% All rights reserved.
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%
%    Redistributions of source code must retain the above copyright notice, this list
%    of conditions and the following disclaimer.
%    Redistributions in binary form must reproduce the above copyright notice, this list 
%    of conditions and the following disclaimer in the documentation and/or other materials 
%    provided with the distribution.
%    Neither the name of the Massachusetts Institute of Technology nor the names of its 
%    contributors may be used to endorse or promote products derived from this software without 
%    specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
% FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
% THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% function d = simple_dijkstra(adj,s)
% 
% n=length(adj);
% d = inf*ones(1,n); % distance s-all nodes
% d(s) = 0;    % s-s distance
% T = 1:n;    % node set with shortest paths not found
% 
% while not(isempty(T))
%     [dmin,ind] = min(d(T));
%     for j=1:length(T)
%         if adj(T(ind),T(j))>0 & d(T(j))>d(T(ind))+adj(T(ind),T(j))
%             d(T(j))=d(T(ind))+adj(T(ind),T(j));
%         end
%     end 
%     T = setdiff(T,T(ind));
%     
% end
