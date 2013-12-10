function [id,od,deg] = degrees_dir(CIJ)
%DEGREES_DIR        Indegree and outdegree
%
%   [id,od,deg] = degrees_dir(CIJ);
%
%   Node degree is the number of links connected to the node. The indegree 
%   is the number of inward links and the outdegree is the number of 
%   outward links.
%
%   Input:      CIJ,    directed (binary/weighted) connection matrix
%
%   Output:     id,     node indegree
%               od,     node outdegree
%               deg,    node degree (indegree + outdegree)
%
%   Notes:  Inputs are assumed to be on the columns of the CIJ matrix.
%           Weight information is discarded.
%
%
%   Olaf Sporns, Indiana University, 2002/2006/2008


% ensure CIJ is binary...
CIJ = double(CIJ~=0);

% compute degrees
id = sum(CIJ,1);    % indegree = column sum of CIJ
od = sum(CIJ,2)';   % outdegree = row sum of CIJ
deg = id+od;        % degree = indegree+outdegree


