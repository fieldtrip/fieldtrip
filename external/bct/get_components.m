function [comps,comp_sizes] = get_components(adj)
%GET_COMPONENTS     connected components
%
%   [comps,comp_sizes] = get_components(adj);
%
%   Returns the components of an undirected graph specified by the binary and 
%   undirected adjacency matrix adj. Components and their constitutent nodes are 
%   assigned the same index and stored in the vector, comps. The vector, comp_sizes,
%   contains the number of nodes beloning to each component.
%
%   Inputs:         adj,    binary and undirected adjacency matrix
%
%   Outputs:      comps,    vector of component assignments for each node
%            comp_sizes,    vector of component sizes
%
%   Note: disconnected nodes will appear as components with a component
%   size of 1
%
%   J Goni, University of Navarra and Indiana University, 2009/2011

%#ok<*ASGLU>

if size(adj,1)~=size(adj,2)
    error('this adjacency matrix is not square');
end

if ~any(adj-triu(adj))
  adj = adj | adj';
end

%if main diagonal of adj do not contain all ones, i.e. autoloops
if sum(diag(adj))~=size(adj,1)
    
    %the main diagonal is set to ones
    adj = adj|speye(size(adj));
end

%Dulmage-Mendelsohn decomposition
[useless,p,useless,r] = dmperm(adj); 

%p indicates a permutation (along rows and columns)
%r is a vector indicating the component boundaries

% List including the number of nodes of each component. ith entry is r(i+1)-r(i)
comp_sizes = diff(r);

% Number of components found.
num_comps = numel(comp_sizes);

% initialization
comps = zeros(1,size(adj,1)); 

% first position of each component is set to one
comps(r(1:num_comps)) = ones(1,num_comps); 

% cumulative sum produces a label for each component (in a consecutive way)
comps = cumsum(comps); 

%re-order component labels according to adj.
comps(p) = comps; 
