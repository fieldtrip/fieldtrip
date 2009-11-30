
function [codebook cdb_labels converg] = lvq1_train(data,labels,codebook_data,codebook_labels,n_iter,alpha,lambda)

% LVQ1_TRAIN  implements the process of LVQ-3 training
% 
% Learning Vector Quantisation-1 (variant LVQ1) - prototype-based
% supervised classification algorithm
%
%  Use as
%    1) [codebook cdb_labels converg] = lvq1_train(data,labels,codebook_data,codebook_labels,n_iter,alpha,lambda)
%    2) As part of LVQ object (see lvq.m)
%
%  INPUT:
%         data   - training data (features)
%         labels - class labels
%         codebook_data,codebook_labels - initial coeebook vectors and their labels
%         n_iter - number of iterations (default: 100 * number of codebook vectors)
%         alpha, lambda - learning rates (see the code and their default settings)
%
%  OUTPUT:
%         codebook, cdb_labels - trained codebook vectors with their labels
%         converg  - log of error convergence
%
%
%  REFERENCES:
%      Somervuo P.; Kohonen T.  Self-Organizing Maps and Learning Vector Quantization for Feature Sequences. 
%                               Neural Processing Letters, Volume 10, Number 2, October 1999, pp.151-159.
%
%  SEE ALSO 
%     lvq3_train.m 
%     lvq.m


%  Pawel Herman, 2009


[n_cod,dim_cod] = size(codebook_data);
[n_data] = size(data,1);

codebook = codebook_data;
converg = sum(sum(codebook.^2)) / (n_cod*dim_cod);

if nargin < 5
    n_iter = 100 * n_cod;
end
if nargin < 6
    alpha = 0.1;
end
if nargin < 7
    lambda = 0.95;
end

for i=1:n_iter    
    old_codebook = codebook;    
    for j=1:n_data
        distances = sum((ones(n_cod,1) * data(j,:) - codebook).^2,2);
        [min_dist,index]  = min(distances);
        if labels(j) == codebook_labels(index)
            codebook(index,:) = update_codebook(codebook(index,:),data(j,:),alpha);
        else
            codebook(index,:) = update_codebook(codebook(index,:),data(j,:),-alpha);
        end     
    end
    alpha = alpha * lambda;    
    converg = [converg sum(sum((codebook - old_codebook).^2)) / (n_cod*dim_cod)];
end

cdb_labels = codebook_labels;

