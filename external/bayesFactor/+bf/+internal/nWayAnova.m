function bf10 = nWayAnova(y,X,varargin)
% Bayes Factor analysis for an N-Way Anova.
% Don't call this directly, use bf.anova instead.
% y = data values
% X = design matrix for ANOVA (indicator vars)  no constant term
%
% Parm/Value pairs
% 'sharedPriors'  - Cell array of vectors indicating which effects (columns
% of X) share the same prior. [{1:nrEffects}]: all effects share the same prior.
% 'options' - Monte Carlo integration and parrallel computation optionss. [bf.options]
% 'scale', - The scale of the prior. A cell array matching the sharedPriors to specify the 
% scale per set of effects that share a prior. Defaults to {sqrt(2)/2}
% 'almostZero' - Improve integration by not including zero itself [0.00001];
% 'continuousIx' - Columns that contain continuous regressors. 
% BK 2018
nrEffects = size(X,2);

p =inputParser;
p.addParameter('sharedPriors',{},@iscell); % Which effects share a prior? A cell array with indices corresponding to columns of X
p.addParameter('options',bf.options);
p.addParameter('scalePerSharedPrior',{sqrt(2)/2},@iscell); 
p.addParameter('continuousIx',[],@isnumeric); % Index (column of X) of continuous covariates
p.addParameter('almostZero',0.00001,@isnumeric); % Integrating from zero can cause problems. Start at something not quite zero. (This is effect size so this is 0 for practical purposes)
p.parse(varargin{:});

if isempty(p.Results.sharedPriors)
    sharedPriors = {1:nrEffects};
else
    sharedPriors = p.Results.sharedPriors;
end
nrDims = numel(sharedPriors);

% Check that the scale matches in each of the dimensions
for i=1:nrDims
    scalesThisDim = size(p.Results.scalePerSharedPrior(i),2);
    if scalesThisDim>1 && scalesThisDim ~=numel(sharedPriors{i}) 
        error('The number of scale elements (%d) does not match the number shared priors (%d)',numel(p.Results.scalePerSharedPrior),nrDims);
    end
end

if nrDims>= p.Results.options.nDimsForMC
    % Use MC Sampling to calculate the integral   
    thisRouderS = @(g)(bf.internal.rouderS(g,y,X,sharedPriors,p.Results.continuousIx,p.Results.options));
    bf10 = bf.internal.mcIntegral(thisRouderS,p.Results.scalePerSharedPrior,p.Results.options);
else
    % Use numeric quadrature integration (potentially more accurate, but
    % slow)
    switch (nrDims)
        case 1
            bf10 = integral(@integrand,p.Results.almostZero,Inf);
        case 2
            bf10 = integral2(@integrand,p.Results.almostZero,Inf,p.Results.almostZero,Inf);
        case 3
            bf10 = integral3(@integrand,p.Results.almostZero,Inf,p.Results.almostZero,Inf,p.Results.almostZero,Inf);
    end
end
    function value = integrand(varargin)
        % Nested function, used by the quadrature integration
        g =cat(1,varargin{:});  % Should be [nrEffects nrValues]  (integratal over dim 2)
        allScales = [p.Results.scalePerSharedPrior{:}];
        pdfG = bf.internal.scaledInverseChiPdf(g,1,reshape(allScales,numel(allScales),1));   
        S = bf.internal.rouderS(g,y,X,sharedPriors,p.Results.continuousIx,p.Results.options);    
        value = S.*prod(pdfG,1);    
    end

end

