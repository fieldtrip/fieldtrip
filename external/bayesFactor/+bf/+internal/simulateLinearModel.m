function out = simulateLinearModel(lm,effectSize,N)
% Simulate data based on a linear model .
% 
% lm = A linearMixedModel that specifies the base model. This function will
% generate simulated data based on the predictions of this model.
% effectSize = 1: use linear model parameters as estimated in the model
%              0: use effect size 0; simulated the null model.
%              other values scale the std dev of the error term. So 2->
%              noise in the simultaed data has half the stdev.
% 
% N = How many observations to generate. Note that N=1 means
% one observation for each of the combinations of all fixed effects.
% Random effects are ignored in the simulation.
% So if lm analyzed a design with 2 groups of subjects, with 4 measuremetns
% each then N=1 would create a table with 8 entries; the minimal complete
% set to do an analysis. N=10 would create a table with 80 entries, all with
% the same mean effect, but noise according to the estimates of the LM.
%
% OUTPUT
% out {1} =table with simulated data.
% out{2} = formula of the linear model.

% Create a complete data table for a single subject who behaves exactly as
% the model (no noise).

errorStd = std(lm.residuals);
% Gett one full set of all combinations of predictor values. 
% Probably not ok for continuous vars, only for categories...
tbl = unique(lm.Variables(:,lm.PredictorNames));
if effectSize==0
    y = zeros(height(tbl),1);
else
    y = predict(lm,tbl);
    errorStd = errorStd/effectSize; % Defining effectsize relative to stdev of the residuals.   
end
tbl.(lm.ResponseName) = y;
tbl = repmat(tbl,[N 1]); % Replicate N times.
% Now add noise as esimated by the model
%res =lm.residuals;
tbl.(lm.ResponseName) = tbl.(lm.ResponseName) + errorStd*randn([height(tbl) 1]);
out = {tbl,char(lm.Formula)};
end
