function S = fastinvc(L)
% L: sparse matrix
% S: covariance matrix, restricted to nonzeros in L
%
% fast 64 bit MEX implementation of the Takahashi equation

[iii,jjj] = find(L);

% diagonal elements
didx = (iii==jjj);
diels = transpose(1:size(L,2));

% set indices for covariance matrix
X = sortrows([[iii(~didx); jjj(~didx); diels] [jjj(~didx); iii(~didx); diels]],[2 1]);

% gives exact result
% S = fastinvex(L,X(:,1)-1,cumsum([0; histc(X(:,2),1:size(L,2))]));
  
% 32 bit stable implementation
% S = fastinvre(L,X(:,1)-1,cumsum([0; histc(X(:,2),1:size(L,2))]));  
 
% 64 bit compatible implementation (preferred)
S = fastinvre64(L,X(:,1)-1,cumsum([0; histc(X(:,2),1:size(L,2))]));  

end