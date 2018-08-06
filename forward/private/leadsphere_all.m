function out=leadsphere_all(xloc,sensorloc,sensorori)
% usage: out=leadsphere_all(xloc,sensorloc,sensorori)

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

[n,nsens]=size(sensorloc); %n=3 m=? 
[n,ndip]=size(xloc);


xlocrep=reshape(repmat(xloc,1,nsens),3,ndip,nsens);
sensorlocrep=reshape(repmat(sensorloc,ndip,1),3,ndip,nsens);
sensororirep=reshape(repmat(sensorori,ndip,1),3,ndip,nsens);

r2=norms(sensorlocrep);
veca=sensorlocrep-xlocrep;
a=norms(veca);
adotr2=dotproduct(veca,sensorlocrep);

gradf1=scal2vec(1./r2.*(a.^2)+adotr2./a+2*a+2*r2);
gradf2=scal2vec(a+2*r2+adotr2./a);
gradf=gradf1.*sensorlocrep-gradf2.*xlocrep;

F=a.*(r2.*a+adotr2);

A1=scal2vec(1./F);
A2=A1.^2;

A3=crossproduct(xlocrep,sensororirep);
A4=scal2vec(dotproduct(gradf,sensororirep));
A5=crossproduct(xlocrep,sensorlocrep);

out=1e-7*(A3.*A1-(A4.*A2).*A5); %%GRB change

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

