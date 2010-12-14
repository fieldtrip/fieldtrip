function  [y] = UGM_Decode_AlphaExpansion(nodePot, edgePot, edgeStruct, decodeFunc,y)
% INPUT
% nodePot(node,class)
% edgePot(class,class,edge) where e is referenced by V,E (must be the same
% between feature engine and inference engine)
%
% OUTPUT
% nodeLabel(node)

[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = edgeStruct.nStates;
maxState = max(nStates);

% Initialize
if nargin < 5
    [junk y] = max(nodePot,[],2);
end
UGM_ConfigurationPotential(y,nodePot,edgePot,edgeStruct.edgeEnds);

% Do Alpha-Expansions until convergence
while 1
    y_old = y;
    
    for s = 1:maxState
        swapPositions = find(y~=s);
        if ~isempty(swapPositions)
            fprintf('Expanding %d\n',s);
            clamped = y;
            clamped(swapPositions) = 0;
            [clampedNP,clampedEP,clampedES] = UGM_makeClampedPotentials(nodePot,edgePot,edgeStruct,clamped);
            
            nClampedNodes = size(clampedNP,1);
            modifiedNP = zeros(nClampedNodes,2);
            for n = 1:nClampedNodes
                modifiedNP(n,:) = [clampedNP(n,s) clampedNP(n,y(swapPositions(n)))];
            end
            nClampedEdges = size(clampedES.edgeEnds,1);
            modifiedEP = zeros(2,2,nClampedEdges);
            for e = 1:nClampedEdges
                n1 = clampedES.edgeEnds(e,1);
                n2 = clampedES.edgeEnds(e,2);
               modifiedEP(:,:,e) = [clampedEP(s,s,e) clampedEP(s,y(swapPositions(n2)),e)
                   clampedEP(y(swapPositions(n1)),s,e) clampedEP(y(swapPositions(n1)),y(swapPositions(n2)),e)];
            end
            clampedES.nStates(:) = 2;
            ytmp = decodeFunc(modifiedNP,modifiedEP,clampedES);
            swapInd = find(ytmp == 1);
            y(swapPositions(swapInd)) = s;
        end
    end
    
    if all(y==y_old)
        break;
    end    
end
