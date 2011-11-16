function Cij=blockwise_conditionalgranger(S,H,Z,cmbindx,n)

for k = 1:size(cmbindx,1)
  % trivariate system
  sel  = cmbindx{k,1};
  nsel = numel(sel);
  siz  = [size(S) 1];
  s123 = reshape(S(sel,:,:), [sqrt(nsel)*[1 1] siz(2:end)]);
  siz  = [size(Z) 1];
  z123 = reshape(Z(sel,:,:), [sqrt(nsel)*[1 1] siz(2:end)]);
  siz  = [size(H) 1];
  h123 = reshape(H(sel,:,:), [sqrt(nsel)*[1 1] siz(2:end)]);
   
  b1  = 1:n{k,1}(1);
  b2  = b1(end)+[1:n{k,1}(2)];
  b3  = b2(end)+[1:n{k,1}(3)];
  b   = [numel(b1) numel(b2) numel(b3)];
  
  % normalization matrix 2->1/3
  tmp1 = -z123(b2,b1)/z123(b1,b1);
  tmp2 = -z123(b3,b1)/z123(b1,b1);
  tmp3 = -(z123(b3,b2)+tmp2*z123(b1,b2))/(z123(b2,b2)+tmp1*z123(b1,b2));
  
  p1  = [eye(b(1)) zeros(b(1),b(2)) zeros(b(1),b(3));
         tmp1        eye(b(2))       zeros(b(2),b(3));
         tmp2      zeros(b(3),b(2))   eye(b(3))    ];
  p2  = [  eye(b(1))       zeros(b(1),b(2)) zeros(b(1),b(3));
         zeros(b(2),b(1))   eye(b(2))       zeros(b(2),b(3));
         zeros(b(3),b(1))   tmp3               eye(b(3))     ];
  
  P    = p2*p1;
  invP = pinv(P);
  
  %bivariate system  
  sel  = cmbindx{k,2};
  nsel = numel(sel);
  siz  = [size(S) 1];
  s2   = reshape(S(sel,:,:), [sqrt(nsel)*[1 1] siz(2:end)]);
  siz  = [size(H) 1];
  h2   = reshape(H(sel,:,:), [sqrt(nsel)*[1 1] siz(2:end)]);
  siz  = [size(Z) 1];
  z2   = reshape(Z(sel,:,:), [sqrt(nsel)*[1 1] siz(2:end)]);
  
  bx1  = 1:n{k,2}(1);
  bx2  = bx1(end)+[1:n{k,2}(2)];
  bx   = [numel(bx1) numel(bx2)];
  
  for kk = 1:size(s2,3)
    HH = h123(:,:,kk)*invP;
    %Q  = [eye(bx(2))               zeros(bx(2),bx(1));
    %      -z2(bx1,bx2)/z2(bx1,bx1)   eye(bx(1))];
    Q  = [eye(bx(1))               zeros(bx(1),bx(2));
          -z2(bx1,bx2)'/z2(bx1,bx1)   eye(bx(2))];
    B  = Q/h2(:,:,kk);
    BB = [B(bx1,bx1)       zeros(b(1),b(2)) B(bx1,bx2);
          zeros(b(2),b(1))   eye(b(2))      zeros(b(2),b(3));
          B(bx2,bx1)       zeros(b(3),b(2)) B(bx2,bx2)];
    FF = BB*HH;
    numer = abs(det(z2(bx1,bx1)));
    denom = abs(det(FF(b1,b1)*z123(b1,b1)*conj(FF(b1,b1))));
    
    Cij(k,kk) = log(numer./denom);
  end
end
