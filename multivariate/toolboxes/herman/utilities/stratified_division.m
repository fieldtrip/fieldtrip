
function [data,design,index] = stratified_division(data,design,ratio)

% STRATIFIED_DIVISION randomly extracts a subset of the input data according to the given ratio.
% The subset is stratified and randomly shuffled (class representatives are appropriately balanced)
%
%  Use as
%      [data,design,index] = stratified_division(data,design,ratio)
%
%  INPUT          
%           data    - input data features (observations in rows)
%           design  - design matrix - class labels assigned to data
%           ratio   - split ratio (between 0 and 1)
%
%  OUTPUT
%           data    - selected data subset 
%           design  - corresponding labels (in the proportion reflecting the
%                                           distribution in the original data set)
%           index   - indices of the extracted subset (w.r.t. original data set)

% Pawel Herman, 2009


labels = unique(design);
data_aux = [];
design_aux = [];
index = [];  

for i=1:length(labels)
    ind = find(design==labels(i));
    ind = ind(randperm(length(ind)));

    index = [index; ind(1:ratio*length(ind))];
end

% shuffle the data
ind = randperm(length(index));
index = index(ind);
data = data(index,:);
design = design(index);

