function K = rbf_ker(X,sig)

%function K = rbf(X,sig)
%
% Computes an rbf kernel matrix from the input data X
%
%INPUTS
% X   =  a matrix containing all samples as rows
% sig = sigma, the kernel width
%
%OUTPUTS
% K = the rbf kernel matrix 
%

 n = size(X,1);
 for i=1:n
     K(i,i)=1;
     for j=1:i-1
         K(i,j)=exp(-norm(X(i,:)-X(j,:))/(2*sig^2));
         K(j,i)=K(i,j);
     end
 end
