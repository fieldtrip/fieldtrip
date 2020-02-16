function X = interaction(Xa,Xb)
% Create all interaction terms from two dummy coded design matrices.
% (See Box II in Rouder et al. 2012)
%
% INPUT
% Xa, Xb = [nrObservations nrA] and [nrObservations nrB]
%           design matrices with matching number of observation (rows)
% OUTPUT
% X = Dummy coded design matrix [nrObservations nrA*nrB]
%
nA = size(Xa,2);
nB = size(Xb,2);
Xa = repmat(Xa,[1 nB]);
Xb = repmat(Xb,[1 nA]);
ix = (repmat(1:nA:nA*nB,[nA 1]) + repmat((0:nA-1)',[1 nB]))';
Xa = Xa(:,ix(:));
X= Xa.*Xb;

end