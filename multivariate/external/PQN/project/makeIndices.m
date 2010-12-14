function [nIndices] = makeIndices(groups)

groups  = groups(:);
nGroups = max(groups);
indices = cell(nGroups,1);
for i=1:nGroups
   indices{i} = find(groups == i);
end
nIndices = zeros(size(indices));
for i=1:nGroups
   nIndices(i) = length(indices{i});
end