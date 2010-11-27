function grad = fif2grad(filename);

% FIF2GRAD constructs a gradiometer definition from a Neuromag *.fif file
% The resulting gradiometer definition can be used by Fieldtrip for forward
% and inverse computations.
%
% Use as
%   grad = fif2grad(filename);
%
% See also READ_HEADER, CTF2GRAD, BTI2GRAD

% Copyright (C) 2004, Joachim Gross
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
% $Id: fif2grad.m 1412 2010-07-15 10:46:19Z vlalit $

% this try-catch construct ensures that missing gradiometer information is
% handeled in a "graceful" way
grad = [];
try
   megmodel('head',[0 0 0],filename);
   [n,s,t]=chaninfo;
   [TY,NA]=chaninfo('type');
   nCoils=sum(TY+1);      % number of coils
   nSensors=length(TY);   % number of sensors
   grad.pnt=zeros(nCoils,3);
   grad.ori=zeros(nCoils,3);
   grad.tra=zeros(nSensors,nCoils);
   grad.unit='cm';
   % define coils
   kCoil=1;
   for k=1:nSensors,
      if (TY(k)==0),  %magnetometer
         grad.pnt(kCoil,:)=100*(t{k}(1:3,4));
         grad.ori(kCoil,:)=t{k}(1:3,3);
         grad.tra(k,kCoil)=1;
         kCoil=kCoil+1;
         grad.label{k}=deblank(s(k,:));
      elseif (TY(k) == 1), %planar gradiometer
         grad.pnt(kCoil,:)=100*(t{k}(1:3,4)-0.008*t{k}(1:3,1)); % multiply with 100 to get cm
         grad.ori(kCoil,:)=t{k}(1:3,3);
         grad.tra(k,kCoil)= -1;
         kCoil=kCoil+1;
         grad.pnt(kCoil,:)=100*(t{k}(1:3,4)+0.008*t{k}(1:3,1));
         grad.ori(kCoil,:)=t{k}(1:3,3);
         grad.tra(k,kCoil)= 1;
         kCoil=kCoil+1;
         grad.label{k}=deblank(s(k,:));
      else
         error('unknown sensor type');
      end
   end
catch
   warning(['gradiometer information could not be extracted from file']);
   warning(['returning an empty grad structure']);
   return;
end
