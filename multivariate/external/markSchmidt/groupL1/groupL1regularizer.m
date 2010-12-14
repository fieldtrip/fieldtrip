function [f] = groupL1regularizer(w,lambda,groups)
f = sum(lambda.*sqrt(accumarray(groups(groups~=0),w(groups~=0).^2)));
end