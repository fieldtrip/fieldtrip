function out=leadsphere_all(xloc,sensorloc,sensorori)
% usage: out=leadsphere_all(xloc,sensorloc,sensorori)

% Copyright (C) 2003, Guido Nolte
% Copyright (C) 2018, Jan-Mathijs Schoffelen
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

% ---Here begins Guido's original implementation
% [n,nsens]=size(sensorloc); %n=3 m=? 
% [n,ndip]=size(xloc);
% 
% 
% xlocrep=reshape(repmat(xloc,1,nsens),3,ndip,nsens);
% sensorlocrep=reshape(repmat(sensorloc,ndip,1),3,ndip,nsens);
% sensororirep=reshape(repmat(sensorori,ndip,1),3,ndip,nsens);
% 
% r2=norms(sensorlocrep);
% veca=sensorlocrep-xlocrep;
% a=norms(veca);
% adotr2=dotproduct(veca,sensorlocrep);
% 
% gradf1=scal2vec(1./r2.*(a.^2)+adotr2./a+2*a+2*r2);
% gradf2=scal2vec(a+2*r2+adotr2./a);
% gradf=gradf1.*sensorlocrep-gradf2.*xlocrep;
% 
% F=a.*(r2.*a+adotr2);
% 
% A1=scal2vec(1./F);
% A2=A1.^2;
% 
% A3=crossproduct(xlocrep,sensororirep);
% A4=scal2vec(dotproduct(gradf,sensororirep));
% A5=crossproduct(xlocrep,sensorlocrep);
% 
% out=1e-7*(A3.*A1-(A4.*A2).*A5); %%GRB change
% 
% return;

% Here begins JMS' adjustment to Guido's implementation, which is much more
% memory friendly with large numbers of dipoles. Essentially, it does not
% work with massively repmatted sensorloc and xloc matrices
[dum,nsens] = size(sensorloc); %n=3 m=?
[dum, ndip] = size(xloc);

r2   = repmat(sqrt(sum(sensorloc.^2)), ndip, 1);
veca = reshape(repmat(sensorloc,ndip,1),3,ndip,nsens)-reshape(repmat(xloc,1,nsens),3,ndip,nsens);
a    = norms(veca);
adotr2 = dotproduct_light(veca,sensorloc');

gradf1 = (1./r2.*(a.^2)+adotr2./a+2*a+2*r2);
gradf2 = (a+2*r2+adotr2./a);

gradf        = zeros(3,ndip,nsens);
gradf(1,:,:) = gradf1*spdiags(sensorloc(1,:)',0,nsens,nsens) - spdiags(xloc(1,:)',0,ndip,ndip)*gradf2;
gradf(2,:,:) = gradf1*spdiags(sensorloc(2,:)',0,nsens,nsens) - spdiags(xloc(2,:)',0,ndip,ndip)*gradf2;
gradf(3,:,:) = gradf1*spdiags(sensorloc(3,:)',0,nsens,nsens) - spdiags(xloc(3,:)',0,ndip,ndip)*gradf2;

A1  = shiftdim(1./(a.*(r2.*a+adotr2)),-1);
A42 = shiftdim(dotproduct_light(gradf,sensorori'),-1).*(A1.^2); % this is the element wise product of A2 and A4 in Guido's original

A3 = crossproduct_light(xloc', sensorori');
A5 = crossproduct_light(xloc',sensorloc');

out        = zeros(size(A3));
out(1,:,:) = 1e-7*(A3(1,:,:).*A1 - A5(1,:,:).*A42);
out(2,:,:) = 1e-7*(A3(2,:,:).*A1 - A5(2,:,:).*A42);
out(3,:,:) = 1e-7*(A3(3,:,:).*A1 - A5(3,:,:).*A42);
return;

function out=crossproduct_light(x,y)
% this should work because the original x is always repmatted along the
% third (original, i.e. sens) dimension, the y along the second (i.e. loc)
m=size(x,1);
k=size(y,1);
out=zeros(3,m,k);
out(1,:,:)=x(:,2)*y(:,3)'-x(:,3)*y(:,2)';
out(2,:,:)=x(:,3)*y(:,1)'-x(:,1)*y(:,3)';
out(3,:,:)=x(:,1)*y(:,2)'-x(:,2)*y(:,1)';
return;

function out=dotproduct_light(x,y)
k = size(y,1);
% this should work because the x is always 3xnlocxnsens, and y is nsensx3
out=shiftdim(x(1,:,:))*spdiags(y(:,1),0,k,k)+shiftdim(x(2,:,:))*spdiags(y(:,2),0,k,k)+shiftdim(x(3,:,:))*spdiags(y(:,3),0,k,k);
return;


function out=crossproduct(x,y)
[n,m,k]=size(x);
out=zeros(3,m,k);
out(1,:,:)=x(2,:,:).*y(3,:,:)-x(3,:,:).*y(2,:,:);
out(2,:,:)=x(3,:,:).*y(1,:,:)-x(1,:,:).*y(3,:,:);
out(3,:,:)=x(1,:,:).*y(2,:,:)-x(2,:,:).*y(1,:,:);
return; 


function out=dotproduct(x,y)
[n,m,k]=size(x);
outb=x(1,:,:).*y(1,:,:)+x(2,:,:).*y(2,:,:)+x(3,:,:).*y(3,:,:);
out=reshape(outb,m,k);
return; 


function result=norms(x)
[n,m,k]=size(x);
resultb=sqrt(x(1,:,:).^2+x(2,:,:).^2+x(3,:,:).^2);
result=reshape(resultb,m,k);
return; 


function result=scal2vec(x)
[m,k]=size(x);
% result=zeros(3,m,k);
% for i=1:3
%     result(i,:,:)=x;
% end
result=reshape(repmat(x(:)', [3 1]), [3 m k]);
return

