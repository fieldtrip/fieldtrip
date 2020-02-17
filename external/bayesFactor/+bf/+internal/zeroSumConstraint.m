function [Xa,Qa] = zeroSumConstraint(X)
% Impose a zero-sum constraint on a dummy coded predictor matrix X
% as in Rouder et al. 2012 . The goal is to make the model estimable
% while equating marginal priors across levels.
%
% By default this is applied to all fixed effects.
%
% INPUT
% X = Dummy coded predictor Matrix [nrObservations nrLevels].
%       By passing a scalar, the function computes only the projection
%       matrix (Qa) for a predictor matrix with that many effects.
% OUTPUT
% Xa = Matrix with zero-sum constraint [nrObservations nrLevels-1]
% Qa -  projection matrix (Xa = X*Qa).
%
% BK - 2018

if isscalar(X)
    % This is the number of effects
    nrEffects = X;
else
    %Design matrix was passed
    nrEffects = size(X,2);
end

%% Follow  Rouder et al. 2012
Sigmaa =eye(nrEffects)- ones([nrEffects nrEffects])/nrEffects;
[eigenVecs,ev]= eig(Sigmaa','vector');
[~,ix] = sort(ev,'desc');
Qa = eigenVecs(:,ix(1:end-1));
%Iaminus1 = eye(nrEffects-1);
%Sigmaa = Qa*Iaminus1*Qa';
if isscalar(X)
    Xa =[];
else
    % Transform design matrix
    Xa = X*Qa;
end

end
