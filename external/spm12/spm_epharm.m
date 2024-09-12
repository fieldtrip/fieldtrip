function [harmonics] = spm_epharm(v,n,a,b,L)
% basis set of prolate spheroidal harmonics  for modelling magnetic
% interfecnce observed by point magnetometers. 
% FORMAT harmonics = spm_epharm(v,n,a,b,L)
%   v               - channel positions (nc x 3 matrix)   
%   n               - channel orientations (nc x 3 matrix)  
%   a               - semi-major axis length (1 x 1 matrix) 
%   b               - semi-minor axis length (1 x 1 matrix)
%   L               - harmonic order (1 x 1 matrix)          
% Output:
%   harmonics       - prolate spheroidal harmonics spanning external space 
%__________________________________________________________________________
% Copyright (C) 2023 Tim Tierney


% focus
c = double(sqrt(a^2-b^2));
T = sum(v.^2,2)+c.^2;
% coordinates
major = sqrt(T + sqrt(T.^2-4*v(:,2).^2*c^2))./sqrt(2);
phi = atan2(v(:,3),v(:,1));

tmp  = v(:,2)./major;

tmp(tmp<-1)=-1;
tmp(tmp>1)=1;
nabla = acos(tmp);

%  cartestian unit vectors to prolate
denom = sqrt(major.^2-c.^2.*cos(nabla).^2);

emajor = major.*sin(nabla).*cos(phi).*n(:,1);
emajor = emajor + major.*sin(nabla).*sin(phi).*n(:,3);
emajor = emajor +sqrt(major.^2-c.^2).*cos(nabla).*n(:,2);
emajor = emajor./denom;

enabla = sqrt(major.^2-c.^2).*cos(nabla).*cos(phi).*n(:,1);
enabla = enabla + sqrt(major.^2-c.^2).*cos(nabla).*sin(phi).*n(:,3);
enabla = enabla -major.*sin(nabla).*n(:,2);
enabla = enabla./denom;

ephi = cos(phi).*n(:,3)-sin(phi).*n(:,1);

% compute metric coefficients
hmajor = sqrt((major.^2-c.^2*cos(nabla).^2)./(major.^2-c.^2));
hnabla = sqrt((major.^2-c.^2*cos(nabla).^2));
hphi = sqrt(major.^2-c.^2).*sin(nabla);


% compute spherical harmonics and derivatives 
[slm,dslmdphi,dslmdnabla] = spm_slm(nabla,phi,L);

% compute prolate harmonic derivatives
pl = plhat(major,a,c,L);
dpl = dplhatdmajor(major,a,c,L);

dpslmdphi = pl.*dslmdphi;
dpslmdnabla = pl.*dslmdnabla;
dpslmdmajor = dpl.*slm;

% sum it all up in a gradient
Gphi = repmat(ephi./hphi,1,size(dpslmdphi,2)).*dpslmdphi;
Gnabla = repmat(enabla./hnabla,1,size(dpslmdphi,2)).*dpslmdnabla;
Gmajor = repmat(emajor./hmajor,1,size(dpslmdphi,2)).*dpslmdmajor;

harmonics = Gmajor+Gphi+Gnabla;
harmonics = bsxfun(@minus,harmonics,mean(harmonics));
harmonics = bsxfun(@rdivide,harmonics,std(harmonics));
end

function [pl] = plhat(major,a,c,li)

pl=zeros(size(major,1),li^2+2*li);
count=1;
for l=1:li
  for  m = -l:l
    pl(:,count) = plmhat(major,a,c,l,m) ;
    count=count+1;
  end
end

end

function [plm] = plmhat(major,a,c,l,m)
c=c+1e-16;
m=abs(m);
lg = gammaln(2*l+1)-gammaln(1)-gammaln(l+1);
lg =lg-gammaln(l-m+1);
num=lg;
denom=lg;

for k = 1:floor(l/2)
  if ((l-2*k-m+1)>0)
    lg = gammaln(2*l-2*k+1)-gammaln(k+1)-gammaln(l-k+1);
    lg =lg-gammaln(l-2*k-m+1);
    lnumtmp = lg+2*k*log(c./major);
    ldenomtmp = lg+2*k*log(c./a);
    num=num+log(1+((-1)^k).*exp(lnumtmp-num));
    denom=denom+log(1+((-1)^k).*exp(ldenomtmp-denom));
  end
end

series = num-denom;
lt1t2 =(l-m)*log(major/a)+ m/2*log((major.^2-c^2)./(a^2-c^2));
plm = exp(lt1t2+series);

end

function [dpl] = dplhatdmajor(major,a,c,li)

dpl=zeros(size(major,1),li^2+2*li);
count=1;
for l=1:li
  for  m = -l:l
    dpl(:,count) = dplmhatdmajor(major,a,c,l,m) ;
    count=count+1;
  end
end

end

function [dplm] = dplmhatdmajor(major,a,c,l,m)
c=c+1e-16;
m=abs(m);

%-u
%--------------------------------------------------------------------------
lt1t2 =(l-m)*log(major/a)+ m/2*log((major.^2-c^2)./(a^2-c^2));
u = exp(lt1t2);

%-dudr
%--------------------------------------------------------------------------
minor = sqrt(major.^2-c^2);
minorref = sqrt(a^2-c^2);
dudr1=log(l-m)+(l-m-1)*log(major)+m*log(minor);
dudr1=dudr1-(l-m)*log(a)-m*log(minorref);
dudr2 =log(m)+(l-m+1)*log(major)+(m-2)*log(minor);
dudr2 = dudr2 -(l-m)*log(a)-m*log(minorref);
dudr = exp(dudr1)+exp(dudr2);

%-v, dvdr
%--------------------------------------------------------------------------
lg = gammaln(2*l+1)-gammaln(1)-gammaln(l+1);
lg =lg-gammaln(l-m+1);
num=lg;
denom=lg;
num2=0;

for k = 1:floor(l/2)
  if ((l-2*k-m+1)>0)
    lg = gammaln(2*l-2*k+1)-gammaln(k+1)-gammaln(l-k+1);
    lg =lg-gammaln(l-2*k-m+1);
    
    lnumtmp = lg+2*k*log(c./major);
    num=num+log(1+((-1)^k).*exp(lnumtmp-num));
    
    ldenomtmp = lg+2*k*log(c./a);
    denom=denom+log(1+((-1)^k).*exp(ldenomtmp-denom));
    
    lnumtmp2 = (-1)^(k)*exp(lg).*(-2*k)*(c./major).^(2*k).*(1./major);
    num2=num2+lnumtmp2;
  end
end

v = exp(num-denom);
dvdr = num2/exp(denom);

%-dplm
%--------------------------------------------------------------------------
dplm = v.*dudr+u.*dvdr;

end