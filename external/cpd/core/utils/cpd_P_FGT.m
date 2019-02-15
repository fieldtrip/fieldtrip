function [P1, Pt1, PX, L]=cpd_P_FGT(X, Y, sigma2, outliers, sigma2_init)

[N,D]=size(X);[M,D]=size(Y);
hsigma=sqrt(2*sigma2);
if outliers==0
   outliers=10*eps; 
end

% FGT parameters
e          = 9;      % Ratio of far field (default e = 10)
K          = round(min([N M 50+sigma2_init/sigma2])); % Number of centers (default K = sqrt(Nx))
p          = 6;      % Order of truncation (default p = 8)

%[e K p]

% computer Pt1 and denomP
[xc , A_k] = fgt_model(Y' , ones(1,M), hsigma, e,K,p);
Kt1 = fgt_predict(X' , xc , A_k , hsigma,e);

ndi=outliers/(1-outliers)*M/N*(2*pi*sigma2)^(0.5*D);
denomP=(Kt1+ndi);
Pt1=1-ndi./denomP;Pt1=Pt1';

% compute P1
[xc , A_k] = fgt_model(X' , 1./denomP, hsigma, e,K,p);
P1 = fgt_predict(Y' , xc , A_k , hsigma,e); P1=P1';

% compute PX
for i=1:D
 [xc , A_k] = fgt_model(X' , X(:,i)'./denomP, hsigma, e,K,p);
 PX(i,:) = fgt_predict(Y' , xc , A_k , hsigma,e); 
end
PX=PX';

L=-sum(log(denomP))+D*N*log(sigma2)/2;