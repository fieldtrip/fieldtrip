function [slope] = getGroupSlope(w,lambda,grad,groups,threshold)

slope = zeros(size(w));

slope(groups==0) = grad(groups==0);

nGroups = max(groups);
for g = 1:nGroups
    groupNorm = norm(w(groups==g));
    gradNorm = norm(grad(groups==g));
    if groupNorm > threshold
       slope(groups==g) = grad(groups==g) + lambda*w(groups==g)/groupNorm; 
    elseif gradNorm > lambda
        slope(groups==g) = grad(groups==g) - lambda*grad(groups==g)/gradNorm;
    end
end