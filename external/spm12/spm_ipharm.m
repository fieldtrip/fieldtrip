function [harmonics] = spm_ipharm(v,n,a,b,L)
% basis set of prolate spheroidal harmonics  for modelling neuronal signals
% observed by point magnetometers. 
% FORMAT harmonics = spm_ipharm(v,n,a,b,L)
%   v               - channel positions (nc x 3 matrix)   
%   n               - channel orientations (nc x 3 matrix)  
%   a               - semi-major axis length (1 x 1 matrix) 
%   b               - semi-minor axis length (1 x 1 matrix)
%   L               - harmonic order (1 x 1 matrix)          
% Output:
%   harmonics       - prolate spheroidal harmonics spanning internal space 
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

%  Cartesian unit vectors to prolate
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

% compute prolate harmonic dervatives
ql = qlhat(major,a,c,L);
dql = dqlhatdmajor(major,a,c,L);

dpslmdphi = ql.*dslmdphi;
dpslmdnabla = ql.*dslmdnabla;
dpslmdmajor = dql.*slm;

% sum it all up in a gradient
Gphi = repmat(ephi./hphi,1,size(dpslmdphi,2)).*dpslmdphi;
Gphi(hphi==0,:)=0;
Gnabla = repmat(enabla./hnabla,1,size(dpslmdphi,2)).*dpslmdnabla;
Gmajor = repmat(emajor./hmajor,1,size(dpslmdphi,2)).*dpslmdmajor;

harmonics = Gmajor+Gphi+Gnabla;
harmonics = bsxfun(@minus,harmonics,mean(harmonics,'omitnan'));
harmonics = bsxfun(@rdivide,harmonics,std(harmonics,'omitnan'));

end

function [ql] = qlhat(major,a,c,li)

ql=zeros(size(major,1),li^2+2*li);
count=1;
for l=1:li
  for  m = -l:l
    ql(:,count) = qlmhat(major,a,c,l,m) ;
    count=count+1;
  end
end

end

function [qlm] = qlmhat(major,a,c,l,m)
m=abs(m);

if c==0
  lt1t2 =(-l-m-1)*log(major/a)+ m/2*log((major.^2-c^2)./(a^2-c^2));
  qlm = exp(lt1t2);
  return
end

k=0;
lg = gammaln((1+l+m)/2+k)+gammaln((2+l+m)/2+k)-gammaln(l+3/2+k);
F1tmp =exp(lg+2*k*log(c./major));
F1 = F1tmp;
F2tmp =exp(lg+2*k*log(c./a));
F2 = F2tmp;

check = [F1tmp;F2tmp];
while(any(check > 1e-32) && k < 5000)
  k=k+1;
  lg = gammaln((1+l+m)/2+k)+gammaln((2+l+m)/2+k)-gammaln(l+3/2+k)-sum(log(1:k));
  F1tmp = exp(lg+2*k*log(c./major));
  F1 =F1+F1tmp;
  F2tmp =exp(lg+2*k*log(c./a));
  F2 = F2+F2tmp;
  check = [F1tmp;F2tmp];
end

lt1t2 =(-l-m-1)*log(major/a)+ m/2*log((major.^2-c^2)./(a^2-c^2));
qlm = exp(lt1t2+log(F1)-log(F2));

end

function [dql] = dqlhatdmajor(major,a,c,li)

dql=zeros(size(major,1),li^2+2*li);
count=1;
for l=1:li
  for  m = -l:l
    dql(:,count) = dqlmhatdmajor(major,a,c,l,m) ;
    count=count+1;
  end
end

end

function [dqlm] = dqlmhatdmajor(major,a,c,l,m)
m=abs(m);
c=c+1e-32;
k=0;
lg = gammaln((1+l+m)/2+k)+gammaln((2+l+m)/2+k)-gammaln(l+3/2+k);
F1tmp =exp(lg+2*k*log(c./major));
F1 = F1tmp;
F2tmp =exp(lg+2*k*log(c./a));
F2 = F2tmp;

check = [F1tmp;F2tmp];
while(any(check > 1e-16) && k < 5000)
  k=k+1;
  lg = gammaln((1+l+m)/2+k)+gammaln((2+l+m)/2+k)-gammaln(l+3/2+k)-sum(log(1:k));
  F1tmp = exp(lg+2*k*log(c./major));
  F1 =F1+F1tmp;
  F2tmp =exp(lg+2*k*log(c./a));
  F2 = F2+F2tmp;
  check = [F1tmp;F2tmp];
end

lt1t2 =(-l-m-1)*log(major/a)+ m/2*log((major.^2-c^2)./(a^2-c^2));
u = exp(lt1t2);
v = F1./F2;


k=0;
lg = gammaln((1+l+m)/2+k)+gammaln((2+l+m)/2+k)-gammaln(l+3/2+k);
g = exp(lg);
dvdmajortmp =g.*(-2*k./major).*((c./major).^(2*k)) ;
dvdmajor = dvdmajortmp;



while(any(dvdmajortmp > 1e-32) && k < 5000)
  k=k+1;
  lg = gammaln((1+l+m)/2+k)+gammaln((2+l+m)/2+k)-gammaln(l+3/2+k)-sum(log(1:k));
  %g = exp(lg);
  hmm = lg+2*k*log(c./major);
  dvdmajortmp =(-2*k./major).*(hmm) ;
  dvdmajor = dvdmajor+dvdmajortmp;
end
dvdmajor = dvdmajor./F2;

minor = sqrt(major.^2-c^2);


lt1 = (-l-m-2)*log(major)+m*log(minor)-(-l-m-1)*log(a)-m/2*log(a^2-c^2);
lt2 = log(m)+(-l-m)*log(major)+(m-2)*log(minor)-(-l-m-1)*log(a)-m/2*log(a^2-c^2);

dudmajor = (-l-m-1)*exp(lt1)+exp(lt2);

dqlm=u.*dvdmajor+v.*dudmajor;

end