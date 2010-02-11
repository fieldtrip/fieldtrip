
% C_bpc - bi-serial point correlation

function [C_bpc,C_wfc,C_ff,C_ff1,C_ff2] = feature_corr(data,design)

% FEATURE_CORR estimates correlation coefficients for a given data set. 
% Feature-feature correlation and bi-serial point correlation feature-class
% is calculated.
% 
%  Use as
%         [C_bpc,C_wfc,C_ff,C_ff1,C_ff2] = feature_corr(data,design) 
%
%  INPUT
%           data    - input data features (observations in rows)
%           design  - design matrix - class labels assigned to data
%
%  OUTPUT
%           C_bpc   - feature-class bi-serial point correlation
%           C_wfc   - weighted linear correlation between features and classes 
%                     (weighting accounts for teh contribution from each class)
%           C_ff    - feature-feature correlation
%           C_ff1/C_ff2 - feature-feature correlation independently in each of two classes

% Pawel Herman, 2009


nfeatures = normalizemeanstd(data);
lab = unique(design);

if length(lab) ~= 2
    error('There are more than two design labeled differently than 1 and 2');
end

ind1 = (design==lab(1));
ind2 = (design==lab(2));
p = length(ind1)/(length(ind1)+length(ind2));
q = length(ind2)/(length(ind1)+length(ind2));

C_bpc = sqrt(p*q) * (mean(nfeatures(ind1,:)) - mean(nfeatures(ind2,:))) ./ std(nfeatures);

%-------------------------------------------------
class1 = ones(length(design),1);
class1(ind2) = 0;
class2 = ones(length(design),1);   
class2(ind1) = 0;
P1 = corr(class1,nfeatures);
P2 = corr(class2,nfeatures);
C_wfc = p * P1(1,:) + q * P2(1,:);
%-------------------------------------------------

C_ff = corr(nfeatures);
C_ff1 = corr(nfeatures(ind1,:));
C_ff2 = corr(nfeatures(ind2,:));

