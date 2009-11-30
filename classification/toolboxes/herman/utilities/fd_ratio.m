
function [fdr diff] = fd_ratio(data,design)

% FD_RATIO calculates Fischer's discriminant ratio (Rayleigh coefficient)
% for two-class data (data1,data2)
% 
%        fdr = abs(mean(data1) - mean(data2)) / (std(data1) + std(data2)) 
%
%  Use as
%         [fdr diff] = fd_ratio(data,design) OR
%         [fdr] = fd_ratio(data,design)
%
%  INPUT
%           data    - input data features (observations in rows)
%           design  - design matrix - class labels assigned to data
%
%  OUTPUT
%           fdr   - Fischer's discriminant ratio
%           diff  - difference between the means (abs(mean(data1)-mean(data2))

% Pawel Herman, 2009

data1 = data(design==1,:);
data2 = data(design==2,:);

m1 = mean(data1,1);
m2 = mean(data2,1);
s1 = std(data1,1);
s2 = std(data2,1);

fdr = abs(m1 - m2)./(s1+s2);
diff = abs(m1 - m2);
