function K = rbf_prim(X,X_test,sig)

n = size(X,1);
m = size(X_test,1);
K = zeros(n,m);

for i=1:n
    for j=1:m
        K(i,j)=exp(-norm(X(i,:)-X_test(j,:))^2/(2*sig^2));
    end
end