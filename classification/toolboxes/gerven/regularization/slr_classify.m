function probs = slr_classify(data,model)
% SLR_CLASSIFY estimates probabilities for given data using multinomial sparse logistic
% regression.
%
% Usage:
%   probs = slr_classify(data,model)
% 
% Input:
%   'data' is a matrix where each row defines the values of the covariates
%   'model' is a cell array that contains the model parameters for each class
%
% Output:
%   'probs' are the class probabilities according to the model
%
% See also: SLR_LEARN
%
% Copyright (c) 2008, Marcel van Gerven
% F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% Remark:
%
% - the bias term is assumed to be part of the data
%
% $Log: slr_classify.m,v $
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%

probs = exp(data * model');

probs = probs ./ repmat(sum(probs,2),1,size(model,1));
