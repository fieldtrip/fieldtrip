function field=meg_forward(dip_par,forwpar)
% calculates the magnetic field of n dipoles
% in a realistic volume conductor
% usage: field=meg_forward(dip_par,forwpar)
%
% input:
% dip_par nx6 matrix where each row contains location (first 3 numbers)
%         and moment (second 3 numbers) of a dipole
% forwpar structure containing all information relevant for this
%         calculation; forwpar is calculated with meg_ini
%         You have here an option to include linear transformations in
%         the forward model by specifying forpwar.lintrafo=A
%         where A is an NxM matrix. Then field -> A field
%         You can use that, e.g., if you can write the forward model
%         with M magnetometer-channels plus a matrix multiplication
%         transforming this to a (eventually higher order) gradiometer.
%
% output:
% field  mxn matrix where the i.th column is the field in m channels
%        of the i.th dipole
%
% note:  No assumptions about units are made (i.e. no scaling factors)
%
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

device_sens=forwpar.device_sens;

field_sens_sphere=getfield_sphere(dip_par,forwpar.device_sens,forwpar.center);
field=field_sens_sphere;clear field_sens_sphere;
if isfield(forwpar,'device_ref')
  field=field-getfield_sphere(dip_par,forwpar.device_ref,forwpar.center);
end
if isfield(forwpar,'device_weights')
  field=field+forwpar.weights*getfield_sphere(dip_par,forwpar.device_weights,forwpar.center);
end

if forwpar.order>0
  coeff=forwpar.coeff_sens;
  if isfield(forwpar,'device_ref')
    coeff=coeff-forwpar.coeff_ref;
  end
  if isfield(forwpar,'device_weights')
    coeff=coeff+forwpar.coeff_weights*forwpar.weights';
  end
  field=field+getfield_corr(dip_par,coeff,forwpar.center,forwpar.order);
end

if isfield(forwpar,'lintrafo');
  field=forwpar.lintrafo*field;
end

return % main function


function field=getfield_sphere(source,device,center)
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
return

function field=getfield_corr(source,coeffs,center,order)
[ndip,ndum]=size(source);
x1=source(:,1:3)-repmat(center',ndip,1);
n1=source(:,4:6);
%correction
if order>0
  scale=10;
  [bas,gradbas]=legs(x1,n1,order,scale);
  nbasis=(order+1)^2-1;
  coeffs=coeffs(1:nbasis,:);
  field=-(gradbas*coeffs)';
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
