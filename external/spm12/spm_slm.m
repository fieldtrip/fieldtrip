function [slm,dslmdphi,dslmdtheta] = spm_slm(theta,phi,li)
% basis set of spehrical harmonics and their (optional) angular derivatives
% as observed by point magnetometers. 
% FORMAT [slm,dslmdphi,dslmdtheta] = spm_slm(theta,phi,li)
%   theta             - colattitude (nchannel x 1 matrix) 
%   phi               - longitude   (nchannel x 1 matrix) 
%   li                - harmonic order (1 x 1 matrix) 
% Output:
%   slm          - spherical harmonics
%   dslmdphi     - derivative with respecte to longitude
%   dslmdtheta   - derivative with respecte to lattitude
%__________________________________________________________________________
% Copyright (C) 2023 Tim Tierney

n=li^2+2*li;
slm=zeros(size(theta,1),n);
dslmdphi=zeros(size(theta,1),n);
dslmdtheta=zeros(size(theta,1),n);

count=1;
for l=1:li
  for  m = -l:l
    a = (-1)^m *sqrt((2*l+1)/(2*pi)* factorial(l-abs(m))/factorial(l+abs(m)));
    % a=1;
    if(m<0)
      L = Plm(cos(theta),l,abs(m));
      slm(:,count) =  a*L.*sin(abs(m)*phi);
      dslmdphi(:,count) =  abs(m)*a*L.*cos(abs(m)*phi);
      L = dPlmdtheta(theta,l,abs(m));
      dslmdtheta(:,count) =  a*L.*sin(abs(m)*phi);
      
    elseif m==0
      L = Plm(cos(theta),l,0);
      slm(:,count) = sqrt((2*l+1)/(4*pi))*L;
      dslmdphi(:,count) = 0;
       L = dPlmdtheta(theta,l,abs(m));
       dslmdtheta(:,count) = sqrt((2*l+1)/(4*pi))*L;
    else
      L = Plm(cos(theta),l,m);
      slm(:,count) =  a*L.*cos(m*phi);
      dslmdphi(:,count) =  (-m)*a*L.*sin(m*phi);
      L = dPlmdtheta(theta,l,m);
      dslmdtheta(:,count) =  a*L.*cos(m*phi);
    end
    count=count+1;
  end
end
end

function [dpl] = dPlmdtheta(theta,l,m)
b= (-1)^m * 2^l ;

dpl =0;
for k = m:l
  tmp = (l+k-1)/2-(0:(l-1));
  val=prod(tmp);
  vals2= prod(l-(0:(k-1)));
  c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
  dpl=dpl+b*c*(m*cos(theta).^(k-m+1).*sin(theta).^(m-1)-(k-m)*sin(theta).^(m+1).*cos(theta).^(k-m-1));
end
dpl(~isfinite(dpl))=0;
end

function [pl] = Plm(x,l,m)
b= (-1)^m * 2^l;
pl =0;
xsq=(1-x.^2).^(m/2);
for k = m:l
  tmp = (l+k-1)/2-(0:(l-1));
  val=prod(tmp);
  vals2= prod(l-(0:(k-1)));
  c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
  pl=pl+b*xsq*c.*x.^(k-m);
end

end

