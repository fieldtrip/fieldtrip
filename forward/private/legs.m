function [basis,gradbasis]=legs(x,dir,n,scale)
% usage: [basis,gradbasis]=legs(x,dir,n,scale)
%
% returns the values and directional derivatives  of (n+1)^2-1 basis functions 
% constructed from spherical harmonics at locations given in x and, for the 
% gradients, for (in general non-normalized) directions given in dir.   
% 
% input: x      set of N locations given as an Nx3 matrix 
%        dir    set of N direction vectors given as an Nx3 matrix 
%                  (dir is not normalized (it hence can be a dipole moment))
%        n       order of spherical harmonics 
%
% output: basis: Nx((n+1)^2-1)  matrix containing in the j.th  row the real 
%                and imaginary parts of r^kY_{kl}(theta,Phi)/(N_{kl}*scale^k) ( (r,theta,phi) 
%                are the spherical coordinates corresponding to  the j.th row in x) 
%                for k=1 to n and l=0 to k 
%                the order is:
%                          real parts for k=1 and l=0,1 (2 terms) 
%                  then    imaginary parts for k=1 and l=1 (1 term) 
%                  then    real parts for k=2 and l=0,1,2 (3 terms) 
%                  then    imaginary parts for k=2 and l=1,2 (2 term) 
%                              etc.
%                   the spherical harmonics are normalized with
%                   N_{kl}=sqrt(4pi (k+l)!/((k-l)!(2k+1)))
%                    the phase does not contain the usual (-1)^l term !!! 
%                   scale is constant preferably set to the avererage radius                   
%
%         gradbasis: Nx((n+1)^2-1) matrix containing in the j.th row the scalar 
%                     product of the gradient of the former with the j.th row of dir

% Copyright (C) 2003, Guido Nolte
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

[n1,n2]=size(x);

comi=sqrt(-1);

normalize=ones(n,n+1);
facto=ones(2*n+2,1);for i=3:2*n+2,facto(i)=facto(i-1)*(i-1);end;
for i=1:n
    for j=1:i+1
        if j==1
            normalize(i,j)=scale^i*i*sqrt(2*facto(i+j-1+1)/facto(i-j+1+1)/(2*i+1));
        else 
            normalize(i,j)=scale^i*i*sqrt(facto(i+j-1+1)/facto(i-j+1+1)/(2*i+1));
        end
    end
end
normalize=reshape(repmat(reshape(normalize,n*(n+1),1),1,n1),n,n+1,n1);


rad=sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
phi=angle(x(:,1)+comi*x(:,2));
costheta=x(:,3)./(rad+eps);

ms=0:n;
ns=1:n;

shiftfactors=zeros(n,n+1);
shiftminusfactors=zeros(n,n+1);
for in=1:n; for im=1:n+1
        shiftfactors(in,im)=in+im-1;
        shiftminusfactors(in,im)=(in+im-2)*(in+im-1);
    end;end;
for in=1:n
    shiftminusfactors(in,1)=1;
end


leg0=zeros(n,n+1,n1);

for i=1:n
 p=legendre(i,costheta);
 leg0(i,1:i+1,:)=p;
end
 leg0=leg0.*reshape(repmat((-1).^ms,n,n1),n,n+1,n1);
ephi=exp(comi*ms'*phi');

leg0=leg0.*reshape(repmat(reshape(ephi,1,(n+1)*n1),n,1),n,n+1,n1)...
         .*reshape(repmat((repmat(rad,1,n).^repmat(ns,n1,1))',n+1,1),n,n+1,n1);

     onesx=zeros(1,n+1,n1);
     onesx(1,1,:)=ones(1,1,n1);
     
legshift=[onesx;leg0(1:n-1,:,:)];
legshiftminus=[-conj(legshift(:,2,:)),legshift(:,1:n,:)].*reshape(repmat(reshape(shiftminusfactors,n*(n+1),1),n1,1),n,n+1,n1);
legshiftplus=-[legshift(:,2:n+1,:),zeros(n,1,n1)];
legshift=legshift.*reshape(repmat(reshape(shiftfactors,n*(n+1),1),n1,1),n,n+1,n1);

dirp=[(dir(:,1)+dir(:,2)/comi)/2,(dir(:,1)-dir(:,2)/comi)/2,dir(:,3)];



gradleg= reshape(repmat(transpose(dirp(:,1)),n*(n+1),1),n,n+1,n1).*legshiftplus...
        +reshape(repmat(transpose(dirp(:,2)),n*(n+1),1),n,n+1,n1).*legshiftminus...
        +reshape(repmat(transpose(dirp(:,3)),n*(n+1),1),n,n+1,n1).*legshift;
    

leg0=leg0./normalize;
gradleg=gradleg./normalize;

    
    basis=zeros(n1,(n+1)^2-1);gradbasis=zeros(n1,(n+1)^2-1);
 for i=1:n
     basis(:,i^2:i^2+i)=(reshape(real(leg0(i,1:i+1,:)),i+1,n1))';
     basis(:,i^2+i+1:(i+1)^2-1)=(reshape(imag(leg0(i,2:i+1,:)),i,n1))';
     gradbasis(:,i^2:i^2+i)=(reshape(real(gradleg(i,1:i+1,:)),i+1,n1))';
     gradbasis(:,i^2+i+1:(i+1)^2-1)=(reshape(imag(gradleg(i,2:i+1,:)),i,n1))';    
 end
 

 
    

return;
