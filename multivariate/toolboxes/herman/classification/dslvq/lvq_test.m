
function labels = lvq_test(data,codebook_data,codebook_labels,weights)

% LVQ_TEST implements the recall phase (test) of LVQ-based classification
% 
%  Use as:
%          1) labels = lvq_test(data,codebook_data,codebook_labels,weights)
%          2) As part of LVQ object (see lvq.m)
%
%  INPUT 
%       data - feature data for testing
%       codebook_data,
%       codebook_labels - codebook vectors and their labels used here for 
%                         classifying new data examples (input data)
%                         (they are the product of learning functions: dslvq, lvq1 or lvq3) 
%       weights - it is possible to assign weighting factors to each feature
%                 component (weights are trained in DSLVQ) - all weights are 1 by default
%
%  OUTPUT
%       labels - class assignments produced for input data
%
% SEE ALSO 
%    dslvq_train.m 
%    lvq1_train.m 
%    lvq3_train.m
%    lvq.m

%  Pawel Herman, 2009

[n_cod,dim_cod] = size(codebook_data);
[n_data] = size(data,1);

if exist('weights','var')
  if isempty(weights), weights = ones(1,dim_cod); end
elseif nargin < 4 
    weights = ones(1,dim_cod);
end

weights(weights < 0) = 0;
weights_gt0_sq = weights.^2;

distances = zeros(n_data,n_cod);
for i=1:n_data
    distances(i,:) = (((ones(n_cod,1) * data(i,:) - codebook_data).^2) * weights_gt0_sq')';
end

[aux ind] = min(distances,[],2);

labels = codebook_labels(ind);