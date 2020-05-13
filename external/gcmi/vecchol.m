function [c] = vecchol(x)
%computes cholesky decomposition of matrix x, using explicit analytic definition if
%size(x,1) < n, otherwise use matlab chol-function

n   = inf;
siz = size(x);
c   = zeros(siz);

siz = size(x);

if numel(siz)==2
    if siz(2)==1
        % 1d response
        c = sqrt(x);
    else
        c = chol(x);
    end
elseif all(siz(2:3)==2)
  %d = x(1,1,:,:).*x(2,2,:,:) - x(1,2,:,:).*x(2,1,:,:);
  A = sqrt(x(:,1,1));
  B = x(:,2,1)./A;
  
  c(:,1,1) = A;
  c(:,1,2) = B;
  c(:,2,2) = sqrt(x(:,2,2)-B.^2);
  
elseif all(siz(2:3)==3),
  A = sqrt(x(:,1,1));
  B = x(:,2,1)./A;
  C = sqrt(x(:,2,2)-B.^2);
  D = x(:,3,1)./A;
  E = (x(:,3,2)-(B.*D))./C;
  F = sqrt(x(:,3,3)-E.^2-D.^2);
  
  c(:,1,1) = A;
  c(:,1,2) = B;
  c(:,1,3) = D;
  c(:,2,2) = C;
  c(:,2,3) = E;
  c(:,3,3) = F;

elseif all(siz(2:3)==4),
  A = sqrt(x(:,1,1));
  B = x(:,2,1)./A;
  C = sqrt(x(:,2,2)-B.^2);
  D = x(:,3,1)./A;
  E = (x(:,3,2)-(B.*D))./C;
  F = sqrt(x(:,3,3)-E.^2-D.^2);
  G = x(:,4,1)./A;
  H = (x(:,4,2)-(B.*G))./C;
  I = (x(:,4,3)-(D.*G)-(E.*H))./F;
  J = sqrt(x(:,4,4)-G.^2-H.^2-I.^2);
  
  c(:,1,1) = A;
  c(:,1,2) = B;
  c(:,1,3) = D;
  c(:,1,4) = G;
  c(:,2,2) = C;
  c(:,2,3) = E;
  c(:,2,4) = G;
  c(:,3,3) = F;
  c(:,3,4) = I;
  c(:,4,4) = J;
  

else
  %write for loop
  for k = 1:siz(1)
    [tmp,dum] = chol(shiftdim(x(k,:,:)));
    if ~dum
      c(k,:,:) = tmp;
    end
  end
end
