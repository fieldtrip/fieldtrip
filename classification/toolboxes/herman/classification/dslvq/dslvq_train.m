
function [codebook, cdb_labels, weights, converg] = dslvq_train(data,labels,codebook_data,codebook_labels,n_iter,alpha,eps,beta,lambda,win_size)

% DSLVQ_TRAIN implements the process of DSLVQ training
% Distinction Sensitive Learning Vector Quantization (DSLVQ)- 
% prototype-based supervised classification algorithm with feature weighting 
%
%  Use as
%    1) [codebook, cdb_labels, weights, converg] = ...
%               dslvq_train(data,labels,codebook_data,codebook_labels,n_iter,alpha,eps,beta,lambda,win_size)
%    2) As part of LVQ object (see lvq.m)
%
%  INPUT:
%         data   - training data (features)
%         labels - class labels
%         codebook_data,codebook_labels - initial coeebook vectors and their labels
%         n_iter - number of iterations (default: 100 * number of codebook vectors)
%         alpha, eps, beta, lambda - learning rates (see the code and their default settings)
%         win_size - size of the update window (default value = 0.3)              
%
%  OUTPUT:
%         codebook, cdb_labels - trained codebook vectors with their labels
%         weights  - weights attached to particular components of feature data (they express their importance)
%         converg  - log of error convergence
%
%
%  REFERENCES:
%      Pregenzer, M; Pfurtscheller, G. Frequency component selection for an EEG-based brain to computer interface. 
%                                      IEEE Transactions on Rehabilitation Engineering. 1999;7(4):413–419.
%
%  SEE ALSO 
%       lvq1_train.m
%       lvq3_train.m 
%       lvq.m


%  Pawel Herman, 2009


[n_cod,dim_cod] = size(codebook_data);
[n_data] = size(data,1);

codebook = codebook_data;
old_codebook = codebook;
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
    beta = 0.1;
end
if nargin < 9
    lambda = 0.95;
end
if nargin < 10
    win_size = 0.3;
end

window = (1-win_size) / (1+win_size);

%Initialize the weight vector
weights	= ones(dim_cod,1);

%while (sum(sum(abs(mu - old_mu))) > 0.1),

for i=1:n_iter
    old_codebook = codebook;

    for j = 1:n_data

        distances = sum((ones(n_cod,1) * data(j,:) - codebook).^2,2);

        % the nearest neighbor classified correctly (positive)
        dp = min(distances(codebook_labels==labels(j)));

        % the nearest one classified incorrectly (negative)
        dn = min(distances(codebook_labels~=labels(j)));

        if isempty(dp) || isempty(dn)
            break
        else
            indp = find(codebook_labels==labels(j) & distances==dp);
            indn = find(codebook_labels~=labels(j) & distances==dn);
            indp = indp(1);
            indn = indn(1);
        end

        distp = abs(data(j,:) - codebook(indp,:));
        distn = abs(data(j,:) - codebook(indn,:));
        nw = (distn-distp)/sum(max(distp,distn));   % nw = (distn-distw)/sum(abs(distp-distn));

        weights = weights + beta * ( nw' - weights);
        weights = weights./sum(abs(weights));   % L1-normalization


        %LVQ3 learning rule for codebooks
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

    %multiplicative decrease
    alpha = lambda * alpha;
    beta  = lambda * beta;
    converg = [converg  sum(sum((codebook - old_codebook).^2)) / (n_cod*dim_cod)];
end

cdb_labels = codebook_labels;
weights = weights'; % row vector