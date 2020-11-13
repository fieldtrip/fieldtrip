function ix = sharedPriorIx(X,terms,sharedPriors)
% Determine which column in the cell array of design matrices X with terms named 'terms'
% should use the same priors as specified by the names of the terms/columns
% in the cell array sharedPriors
% X - cell array of design matrices
% terms - cell array with names of each term/column in X
% sharedPriors - cell array of names/terms that should share priors.
% OUTPUT
% ix = columns of the concatenated design matrix [X{:}] that should share
% their priors.
% 
% Used by bf.anova.m
%

% Based on the names in sharedPriors, find indices of columns that need to
% share priors.
% Assign groups of effects to use the same prior on effect size.
ix = cell(1,numel(sharedPriors));
[ix{:}] = deal([]);
soFar  =0;
if isempty(sharedPriors)
    ix = {};
else
    for i=1:numel(terms)
        match = cellfun(@any,cellfun(@(x)(strcmp(terms{i},x)),sharedPriors,'UniformOutput',false));
        if ~any(match)
            error(['Shared priors not defined for ' terms{i}]);
        end
        nrInThisTerm  = size(X{i},2);
        ix{match} = cat(2,ix{match},soFar+(1:nrInThisTerm));
        soFar = soFar+nrInThisTerm;
    end
end
end