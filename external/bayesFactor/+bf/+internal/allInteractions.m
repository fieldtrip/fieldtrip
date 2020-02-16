function interactions = allInteractions(factorNames)
% Given a cell array of factor names, create all pairwise combinations
% of factors to represent interactions
% INPUT
% factorNames  - cell array of names
% OUTPUT
% interactionNames - cell array of interaction names
%
% allInteractions({'a','b'}) -> {'a:b'}
%
% BK  -Nov 2018

cntr=0;
nrFactors = numel(factorNames);
nrInteractions = nrFactors*(nrFactors-1)-1;
interactions =cell(1,nrInteractions);
for i=1:nrFactors
    for j=(i+1):nrFactors
        cntr= cntr+1;
        interactions{cntr} = [factorNames{i} ':' factorNames{j}];
    end
end
end
