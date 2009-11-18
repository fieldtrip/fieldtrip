function ip = anova(data,design)
% ANOVA computes 1-pvalues that data in different classes belong to the
% same group. I.e., the closer to one, the better.
%
%   ip = anova(data,design)
%
%   ip = 1 - pvalue
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: anova.m,v $
%
    ip = zeros(1,size(data,2));
    for j=1:size(data,2)        
        ip(j) = 1 - anova1(data(:,j),design(:,1),'off');
    end
    
end
