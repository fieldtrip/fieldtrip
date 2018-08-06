function gt = gtom(adj,numSteps)
%GTOM       Generalized topological overlap measure
%
%   gt = gtom(adj,numSteps);
%
%   The m-th step generalized topological overlap measure (GTOM) quantifies
%   the extent to which a pair of nodes have similar m-th step neighbors.
%   Mth-step neighbors are nodes that are reachable by a path of at most
%   length m.
%
%   This function computes the the M x M generalized topological overlap
%   measure (GTOM) matrix for number of steps, numSteps. 
%
%   Inputs:       adj,    adjacency matrix (binary,undirected)
%            numSteps,    number of steps
%
%   Outputs:       gt,    GTOM matrix
%
%   NOTE: When numSteps is equal to 1, GTOM is identical to the topological
%   overlap measure (TOM) from reference [2]. In that case the 'gt' matrix
%   records, for each pair of nodes, the fraction of neighbors the two
%   nodes share in common, where "neighbors" are one step removed. As
%   'numSteps' is increased, neighbors that are furter out are considered.
%   Elements of 'gt' are bounded between 0 and 1.  The 'gt' matrix can be
%   converted from a similarity to a distance matrix by taking 1-gt.
%
%   References: [1] Yip & Horvath (2007) BMC Bioinformatics 2007, 8:22
%               [2] Ravasz et al (2002) Science 297 (5586), 1551.
%
%   J Goni, University of Navarra and Indiana University, 2009/2011

%#ok<*ASGLU>

%initial state for bm matrix;
bm = adj;
bmAux = bm;
numNodes = size(adj,1);

if (numSteps > numNodes)
    disp('warning, reached maximum value for numSteps. numSteps reduced to adj-size')
    numSteps = numNodes;
end

if (numSteps == 0)
    %GTOM0
    gt = adj;
else
    
    for steps = 2:numSteps
        for i = 1:numNodes
            
            %neighbours of node i
            [neighRow,neighColumn] = find(bm(i,:)==1); 
            
            %neighbours of neighbours of node i
            [neighNeighRow,neighNeighColumn] = find(bm(neighColumn,:)==1);
            newNeigh = setdiff(unique(neighNeighColumn),i);
            
            %neighbours of neighbours of node i become considered node i neighbours
            bmAux(i,newNeigh) = 1;
            
            %keep symmetry of matrix
            bmAux(newNeigh,i) = 1;
        end
        %bm is updated with new step all at once
        bm = bmAux;
        
    end
    
    clear bmAux newNeigh;
    
    %numerators of GTOM formula
    numeratorMatrix = bm*bm + adj + speye(numNodes,numNodes);
    
    %vector containing degree of each node
    bmSum=sum(bm);  
    clear bm;
    
    denominatorMatrix = -adj + min(repmat(bmSum,numNodes,1),repmat(bmSum',1,numNodes)) + 1;
    gt = numeratorMatrix ./ denominatorMatrix;
end
