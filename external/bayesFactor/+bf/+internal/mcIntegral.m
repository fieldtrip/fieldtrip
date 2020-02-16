function v = mcIntegral(fun,scale,options)
% Monte Carlo integration
%
% INPUT
% fun -  The function to integrate. This should be specified as a
%       function_handle that takes a single input (g)
% scale - a cell array of scales for each of the prior effect size dimensions 
%
% options - A struct with options.  (see bf.options)
% OUTPUT
% v -  The value of the integral. (Typically the BF10).


%% Sample effect sizes (g)  for each of the dimensions (i.e. effects
%  that share a prior effect size) according to the prior probability of
% those g (using the scaled inverse chi)
gRange =  (options.minG:options.stepG:options.maxG);
nrDims = numel(scale);
g =  nan(nrDims,options.nrSamples);
for i=1:nrDims
    pdf   = bf.internal.scaledInverseChiPdf(gRange,1,scale{i}); % Weights for randsample
    g(i,:) = randsample(gRange,options.nrSamples,true,pdf);
end

%% Evaluate the function at these g values
bf10Samples = fun(g);
v = mean(bf10Samples); % Expectation value ~ integral.
end
