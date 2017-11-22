function dT = ft_deriv(T)

% derivative along the columns, where the 2:end-1 elements are the average
% of the n-1 and n+1 differences

[i1,i2]   = size(T);
dT(i1,i2) = 0;

dT(1,:)  = T(2,:)  - T(1,:);
dT(i1,:) = T(i1,:) - T(i1-1,:);
dT(2:(i1-1),:) = (T(3:i1,:)-T(1:(i1-2),:))./2;