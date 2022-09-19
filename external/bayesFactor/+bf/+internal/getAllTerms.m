function allTerms = getAllTerms(lm)
% Returns a cell array with the names of all terms in the linear model.
allTerms = lm.CoefficientNames;
allTerms(strcmpi(allTerms,'(Intercept)')) = []; % Intercept not neeed
% In the linearmixedmodel the factors have _level appended in their
% names. Even without grouping. I remove this here.
allTerms= regexprep(allTerms,'_(?<level>[\w\d]+)\>','');
allTerms = unique(allTerms);
end
