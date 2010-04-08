
function [codebook cdb_labels converg] = lvq3_train(data,labels,codebook_data,codebook_labels,n_iter,alpha,eps,lambda,win_size)

% LVQ3_TRAIN implements the process of LVQ-3 training
%
% Learning Vector Quantisation-3 (variant LVQ3) - prototype-based
% supervised classification algorithm
%
%  Use as
%    1) [codebook cdb_labels converg] = ...
%                       lvq3_train(data,labels,codebook_data,codebook_labels,n_iter,alpha,eps,lambda,win_size)
%    2) As part of LVQ object (see lvq.m)
%
%  INPUT:
%         data   - training data (features)
%         labels - class labels
%         codebook_data,codebook_labels - initial coeebook vectors and their labels
%         n_iter - number of iterations (default: 100 * number of codebook vectors)
%         alpha, eps, lambda - learning rates (see the code and their default settings)
%         win_size - size of the update window (default value = 0.3)              
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
%    lvq1_train.m 
%    lvq.m


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
    eps = 0.3;
end
if nargin < 8
    lambda = 0.95;
end
if nargin < 9
    win_size = 0.3;
end

window = (1-win_size) / (1+win_size);

for i=1:n_iter
    old_codebook = codebook;    
    for j=1:n_data
        distances = sum((ones(n_cod,1) * data(j,:) - codebook).^2,2);        
        [min_dist1,cod1]  = min(distances);
        distances(cod1) = Inf;
        [min_dist2,cod2] = min(distances);
        
        if min_dist1/min_dist2 > window

            if labels(j)==codebook_labels(cod1)
                if labels(j)==codebook_labels(cod2)
                    codebook(cod1,:) = update_codebook(codebook(cod1),data(j,:),alpha * eps);
                    codebook(cod2,:) = update_codebook(codebook(cod1),data(j,:),alpha * eps);
                else
                    codebook(cod1,:) = update_codebook(codebook(cod1),data(j,:),alpha);
                    codebook(cod2,:) = update_codebook(codebook(cod2),data(j,:),-alpha);
                end
            elseif labels(j)==codebook_labels(cod2)
                codebook(cod1,:) = update_codebook(codebook(cod1),data(j,:),-alpha);
                codebook(cod2,:) = update_codebook(codebook(cod2),data(j,:),alpha);
            end

        end
    end
    alpha = alpha * lambda;
    converg = [converg sum(sum((codebook - old_codebook).^2)) / (n_cod*dim_cod)];
end

cdb_labels = codebook_labels;