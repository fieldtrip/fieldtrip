function [w] = groupSoftThreshold(w,alpha,lambda,groups)
    nGroups = max(groups);
    reg = sqrt(accumarray(groups(groups~=0),w(groups~=0).^2));
    for g = 1:nGroups
        if reg(g) ~= 0
            w(groups==g) = (w(groups==g)/reg(g))*max(0,reg(g)-lambda(g)*alpha);
        end
    end
end