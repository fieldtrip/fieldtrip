function dif = meandiff(data,design)
% MEANDIFF computes the l2 norm of the difference between the means of the data of
% different classes w.r.t. the global mean per feature
%
%   dif = meandiff(data,design)
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: meandif.m,v $
%

    nclasses = max(design(:,1));
    nfeatures = size(data,2);

    dif = zeros(nclasses,nfeatures);
    for j=1:nclasses        
        dif(j,:) = mynanmean(data(design==j,:));
    end
    dif = dif - repmat(mynanmean(data),[nclasses 1]);
    dif = sqrt(sum(dif.^2));

end
