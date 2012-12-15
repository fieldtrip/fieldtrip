function loads = sb_calc_ven_loads(pos,dir,ven_nd,node);

% SB_CALC_VEN_LOADS
%
% $Id$

%aref setzen
aref = 20;
%lambda setzen
lambda = 10e-6;
%r setzen
r = 1;

loads = zeros(size(pos,1),size(ven_nd,2));
for i=1:size(pos,1);
    x = bsxfun(@minus,node(ven_nd(i,ven_nd(i,:)~=0),:),pos(i,:))./aref;
    X = zeros(9,size(x,1));
    X(1:3:7,:) = ones(3,size(x,1));
    X(2:3:8,:) = x';
    X(3:3:9,:) = (x.^2)';
    T = zeros(9,1);
    T(1:3:7) = 0;
    T(2:3:8) = dir(i,:)./aref;
    T(3:3:9) = 0;
    W = zeros(size(x,1)*3,size(x,1));
    W(1:size(x,1),:) = diag(x(:,1).^r);
    W(size(x,1)+1:2*size(x,1),:) = diag(x(:,2).^r);
    W(size(x,1)*2+1:3*size(x,1),:) = diag(x(:,3).^r);
    tmp = ((X')*X + lambda*(W')*W) \ X'*T;
    loads(i,ven_nd(i,:)~=0) = tmp;
end
end
