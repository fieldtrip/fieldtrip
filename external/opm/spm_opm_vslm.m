function [vSlm] = spm_opm_vslm(S)
% cartesian real regular/irregular vector spherical harmonics
% FORMAT vSlm   =  spm_opm_vslm(S)
%   S               - input structure
% Optional fields of S:
% SENSOR LEVEL INFO
%   S.D         - SPM MEEG object      - Default: specify S.D or, S.v and S.o
%   S.v         - optional positions               - Default: same as S.D
%   S.o         - optional orientations            - Default: same as S.D
%   S.or        - optional origin offset           - Default = [0,0,0]
%   S.reg       - regular or irregular (boolean)   - Default: 1
%   S.scale     - scale harmonic for stabilty      - Default: 1
%   S.li        - order of harmonic                - Default: 1
% Output:
%  vSlm            - matrix of vector spherical harmonic (n x (li^2+2*l))
%__________________________________________________________________________
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging
%
% Tim Tierney
% $Id$

% sp = spm_mesh_sphere(5);
% 
% 
% S=[]
% S.li=3;
% S.reg=1;
% S.v=sp.vertices*50;
% S.o=sp.vertices;
% H2= spm_opm_vslm(S);
%
%  for i =1:size(H2,2)
%  f=figure()
%  p= [];
%  p.vertices=sp.vertices;
%  p.faces=sp.faces;
%  p.EdgeColor='none';
%  col= H2(:,i);
%  patch(p,'FaceVertexCData',col,'FaceColor','interp');
%  view([59,28])
%  end
%- handle arguments
%----------------------------------------------------------------------
if ~isfield(S, 'li')
S.li=1;
end 
%- Order 
%----------------------------------------------------------------------

%- MEG object or matrices
%----------------------------------------------------------------------

if ~isfield(S, 'D')
    x=S.v(:,1);
    y=S.v(:,2);
    z=S.v(:,3);
    nx = S.o(:,1);
    ny = S.o(:,2);
    nz = S.o(:,3);
    v=S.v;
else
    p = sensors(S.D,'MEG');
    o=[p.chanori];
    
    nx = o(:,1);
    ny = o(:,2);
    nz = o(:,3);
    
    x = p.chanpos(:,1);
    y = p.chanpos(:,2);
    z = p.chanpos(:,3);
    
end

%- set the origin offset
%----------------------------------------------------------------------
if ~isfield(S, 'or')
    S.or = [0,0,0];
end

x = (x)-S.or(1);
y = (y)-S.or(2);
z = (z)-S.or(3);
v = [x,y,z];

%- regular or irregular harmonics
%----------------------------------------------------------------------
if ~isfield(S, 'reg')
    S.reg=1;
end
reg=S.reg;

%- scaling of coordinate system to prevent infinities
%----------------------------------------------------------------------

if ~isfield(S, 'scale')
    S.scale=1;
end

if(S.scale)
    % need high precision to prevent overflow
    x= double(x);
    y= double(y);
    z= double(z);
    
    sc= max(sqrt(x.^2+y.^2+z.^2));
    x= x/sc;
    y=y/sc;
    z=z/sc;
    v = [x,y,z];
end
%- prepare 
%----------------------------------------------------------------------

n=S.li^2+2*S.li;
r = sqrt(x.^2+y.^2+z.^2);
rbar = mean(r);


xy = x./(x.^2+y.^2);
yx = y./(x.^2+y.^2);
atyx = atan2(y,x);

Slm = slm(v,S.li);
Slmdx=zeros(size(x,1),n);
Slmdy=zeros(size(x,1),n);
Slmdz=zeros(size(x,1),n);
vSlm=zeros(size(x,1),n);
%- dSdX dSdY dSdZ
%----------------------------------------------------------------------
count=1;
for l=1:S.li
    for m=-l:l
        
        a = (-1)^m * sqrt((2*l+1)/(2*pi)* factorial(l-abs(m))/factorial(l+abs(m)));
        u = m*atyx;
        um =abs(m)*atyx;
        L = plm(z./r,l,abs(m));
        [Xlm,Ylm,Zlm]= dplm(v,l,abs(m));
       %[sum(isinf(Xlm)|isnan(Xlm)),sum(isinf(Ylm)|isnan(Ylm)),sum(isinf(Zlm)|isnan(Zlm)),count]
       
          Xlm(isinf(Xlm))=0;
          Ylm(isinf(Ylm))=0;
          Zlm(isinf(Zlm))=0;
          
          
        
        if(m<0)
            %z
            t1= a*sin(um).*Zlm;
            t1(isnan(t1))=0;   
            %y
            t2= a*L*abs(m).*cos(um).*xy+a*sin(um).*Ylm;
            t2(isnan(t2))=0;
            %x
            t3= -a*L*abs(m).*cos(um).*yx+a*sin(um).*Xlm;
            t3(isnan(t3))=0;
        end
  
        if(m==0)
            %z
            t1= sqrt((2*l+1)/(4*pi))*Zlm;
            t1(isnan(t1))=0;  
            %y
            t2= sqrt((2*l+1)/(4*pi))*Ylm; 
            t2(isnan(t2))=0; 
            %x
            t3= sqrt((2*l+1)/(4*pi))*Xlm;
            t3(isnan(t3))=0;
        end
        
        if(m>0)
            %z
            t1= a*cos(u).*Zlm;
            t1(isnan(t1))=0;
            %y
            t2= -a*L*m.*sin(u).*xy+a*cos(u).*Ylm;
            t2(isnan(t2))=0;         
            %x
            t3=a*L*m.*sin(u).*yx+a*cos(u).*Xlm;
            t3(isnan(t3))=0;
        end
  
        if(reg)
            Slmdz(:,count) = t1.*r.^(l)+l*z.*Slm(:,count).*r.^(l-2);
            Slmdy(:,count) = t2.*r.^(l)+l*y.*Slm(:,count).*r.^(l-2);
            Slmdx(:,count) = t3.*r.^(l)+l*x.*Slm(:,count).*r.^(l-2);     
        else
            Slmdz(:,count) = t1./r.^(l+1)-(l+1)*z.*Slm(:,count)./r.^(l+3);
            Slmdy(:,count) = t2./r.^(l+1)-(l+1)*y.*Slm(:,count)./r.^(l+3);
            Slmdx(:,count) = t3./r.^(l+1)-(l+1)*x.*Slm(:,count)./r.^(l+3);
        end
        count=count+1;
    end
end

%- cleanup
%----------------------------------------------------------------------
for i = 1:n
     
vSlm(:,i)=Slmdz(:,i).*nz+Slmdy(:,i).*ny+Slmdx(:,i).*nx;
%[[sum(isnan(Slmdx(:,i))), sum(isnan(Slmdy(:,i))), sum(isnan(Slmdz(:,i))) sum(isnan(vSlm(:,i))) i];
%[sum(isinf(Slmdx(:,i))), sum(isinf(Slmdy(:,i))), sum(isinf(Slmdz(:,i))) sum(isinf(vSlm(:,i))) i]]

end


end
function [Xlm,Ylm,Zlm] = dplm(v,l,m)
x=v(:,1);
y=v(:,2);
z=v(:,3);
r= sqrt(x.^2+y.^2+z.^2);
b= (-1)^m * 2^l;
Xlm=0;
Ylm=0;
Zlm=0;
for k = m:l
    val=prod((l+k-1)/2-(0:(l-1)));
    vals2= prod(l-(0:(k-1)));
    c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
    
    numx = -x.*z.^(k-m).*(k-m).*(x.^2+y.^2).^(m/2) + (m*x.*z.^(k-m+2)).*(x.^2+y.^2).^((m-2)/2);
    numx(isinf(numx))=0;
    numx(isnan(numx))=0;
    
    Xlm=Xlm+b*c*(numx)./(r.^(2+k));
    
    numy = -y.*z.^(k-m).*(k-m).*(x.^2+y.^2).^(m/2) + (m.*y.*z.^(k-m+2)).*(x.^2+y.^2).^((m-2)/2);
    numy(isinf(numy))=0;
    numy(isnan(numy))=0;
    
    Ylm=Ylm+b*c*(numy)./(r.^(2+k));
    
    numz= z.^(k-m-1).*(k-m).*(x.^2+y.^2).^((m+2)/2) + (-m*z.^(k-m+1)).*(x.^2+y.^2).^(m/2);
    numz(isinf(numz))=0;
    numz(isnan(numz))=0;
    
    Zlm=Zlm+b*c*(numz)./(r.^(2+k));     
end
end

function [Slm] = slm(v,li)
x=v(:,1);
y=v(:,2);
z=v(:,3);
r= sqrt(x.^2+y.^2+z.^2);
n=li^2+2*li;
Slm=zeros(size(v,1),n);
count=1;
for l=1:li
    for  m = -l:l
        a = (-1)^m *sqrt((2*l+1)/(2*pi)* factorial(l-abs(m))/factorial(l+abs(m)));
       % a=1;
        if(m<0)
            L = plm(z./r,l,abs(m));
            Slm(:,count) =  a*L.*sin(abs(m)*atan2(y,x));
        elseif m==0
            L = plm(z./r,l,0);
            Slm(:,count) = sqrt((2*l+1)/(4*pi))*L;
        else
            L = plm(z./r,l,m);
            Slm(:,count) =  a*L.*cos(m*atan2(y,x));

        end
        count=count+1;
    end
end

    
end
function [pl] = plm(x,l,m)
b= (-1)^m * 2^l;
pl =0;
for k = m:l
    val=prod((l+k-1)/2-(0:(l-1)));
    vals2= prod(l-(0:(k-1)));
    c = (factorial(k)/factorial(k-m)) * vals2/factorial(k) * val/factorial(l);
    pl=pl+(b*(1-x.^2).^(m/2))*c.*x.^(k-m);
end
end