function forwpar=meg_ini(vc,center,order,sens,refs,gradlocs,weights)
% initializes MEG-forward calculation
% usage: forwpar=meg_ini(vc,center,order,sens,refs,gradlocs,weights)
%
% input:
% vc:  Nx6 matrix; N is the number of surface points
%      the first three numbers in each row  are the location
%      and the second three are the orientation of the surface
%      normal
% center: 3x1 vector denoting the center of volume the conductor
% order: desired order of spherical spherical harmonics;
%        for 'real' realistic volume conductors order=10 is o.k
% sens:  Mx6 matrix containing sensor location and orientation,
%        format as for vc
% refs: optional argument.  If provided, refs contains the location and oriantion
%       (format as sens) of additional sensors which are subtracted from the original
%        ones. This makes a gradiometer. One can also do this with the
%        magnetometer version of this program und do the subtraction outside this program,
%        but the gradiometer version is faster.
% gradlocs, weights: optional two arguments (they must come together!).
%                    gradlocs are the location of additional  channels (e.g. to calculate
%                    a higher order gradiometer) and weights. The i.th row in weights contains
%                    the weights to  correct if the i.th cannel. These extra fields are added!
%                    (has historical reasons).
%
% output:
% forpwar: structure containing all parameters needed for forward calculation
%
% note: it is assumed that locations are in cm.

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

if nargin==4
  if order>0;
    coeff_sens=getcoeffs(sens,vc,center,order);
    forwpar=struct('device_sens',sens,'coeff_sens',coeff_sens,'center',center,'order',order);
  else
    forwpar=struct('device_sens',sens,'center',center,'order',order);
  end
elseif nargin==5
  if order>0;
    coeff_sens=getcoeffs(sens,vc,center,order);
    coeff_refs=getcoeffs(refs,vc,center,order);
    forwpar=struct('device_sens',sens,'device_ref',refs,'coeff_sens',coeff_sens,'coeff_ref',coeff_refs,'center',center,'order',order);
  else
    forwpar=struct('device_sens',sens,'device_ref',refs,'center',center,'order',order);
  end
elseif nargin==7;
  if order>0;
    coeff_sens=getcoeffs(sens,vc,center,order);
    coeff_refs=getcoeffs(refs,vc,center,order);
    coeff_weights=getcoeffs(gradlocs,vc,center,order);
    forwpar=struct('device_sens',sens,'device_ref',refs,'coeff_sens',coeff_sens,'coeff_ref',coeff_refs,'center',center,'order',order,'device_weights',gradlocs,'coeff_weights',coeff_weights,'weights',weights);
  else
    forwpar=struct('device_sens',sens,'device_ref',refs,'center',center,'order',order,'device_weights',gradlocs,'weights',weights);
  end
  
else
  error('you must provide 4,5 or 7 arguments');
end

return % main function


function coeffs=getcoeffs(device,vc,center,order)
[ndip,ndum]=size(vc);
[nchan,ndum]=size(device);
x1=vc(:,1:3)-repmat(center',ndip,1);
n1=vc(:,4:6);
x2=device(:,1:3)-repmat(center',nchan,1);
n2=device(:,4:6);
scale=10;
nbasis=(order+1)^2-1;
[bas,gradbas]=legs(x1,n1,order,scale);
bt=leadsphere_all(x1',x2',n2');
n1rep=reshape(repmat(n1',1,nchan),3,ndip,nchan);
b=dotproduct(n1rep,bt);
ctc=gradbas'*gradbas;
warning('OFF', 'MATLAB:nearlySingularMatrix');
coeffs=inv(ctc)*gradbas'*b;
warning('ON', 'MATLAB:nearlySingularMatrix');
return

function field=getfield(source,device,coeffs,center,order)
[ndip,ndum]=size(source);
[nchan,ndum]=size(device);
x1=source(:,1:3)-repmat(center',ndip,1);
n1=source(:,4:6);
x2=device(:,1:3)-repmat(center',nchan,1);
n2=device(:,4:6);
%spherical
bt=leadsphere_all(x1',x2',n2');
n1rep=reshape(repmat(n1',1,nchan),3,ndip,nchan);
b=dotproduct(n1rep,bt);
field=b';
%correction
if order>0
  scale=10;
  [bas,gradbas]=legs(x1,n1,order,scale);
  nbasis=(order+1)^2-1;
  coeffs=coeffs(1:nbasis,:);
  fcorr=gradbas*coeffs;
  field=field-fcorr';
end
return

function out=crossproduct(x,y)
% usage: out=testprog(x,y)
% testprog calculates the cross-product of vector x and y
[n,m,k]=size(x);
out=zeros(3,m,k);
out(1,:,:)=x(2,:,:).*y(3,:,:)-x(3,:,:).*y(2,:,:);
out(2,:,:)=x(3,:,:).*y(1,:,:)-x(1,:,:).*y(3,:,:);
out(3,:,:)=x(1,:,:).*y(2,:,:)-x(2,:,:).*y(1,:,:);
return

function out=dotproduct(x,y)
% usage: out=dotproduct(x,y)
% testprog calculates the dotproduct of vector x and y
[n,m,k]=size(x);
outb=x(1,:,:).*y(1,:,:)+x(2,:,:).*y(2,:,:)+x(3,:,:).*y(3,:,:);
out=reshape(outb,m,k);
return

function result=norms(x)
[n,m,k]=size(x);
resultb=sqrt(x(1,:,:).^2+x(2,:,:).^2+x(3,:,:).^2);
result=reshape(resultb,m,k);
return
