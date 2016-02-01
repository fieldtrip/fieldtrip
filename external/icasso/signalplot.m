function h_line=signalplot(S)
%function h_line=signalplot(S)
%
%PURPOSE
%
%To plot several signals, e.g., estimated sources into the same
%graphic axis; a simplified alternative to MATLAB's native subplot
%command.  
%
%EXAMPLE OF BASIC USE
%
%Get the "best" estimated sources from Icasso result struct sR and
%plot them:
%
% [iq,A,W,S]=icassoResult(sR);
% signalplot(S);
%
%INPUTS
%
% S any (MxN matrix) of M signals   
%
%OUTPUT
%
% h_line (vector) graphic handles to the line objects
%
%DETAILS
%
%Plots signals into the active graphic axis. 
%
%If there are more than one signal the function rescales all signals
%between 0 and 1, and sets y-axis tick marks to 1,2,...,M. 
%
%If there is only one signal, the function just plots it using MATLAB plot
%function.  
%
%SEE ALSO
% icassoShow

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.2 100105 johan

[M N]=size(S);

% Default case is to put signals in y-coordinates 1,2,...,M
ycoord=[1:M]';

GAP=0.03;

cla reset;
if M>1,
  % Set x-coordinates
  x=repmat(1:N,M,1);
  %Compute scaling parameters
  mn=repmat(min(S,[],2),1,N);  
  mx=repmat(max(S,[],2),1,N); e=mx-mn;
  e(e==0)=1;
  S=(S-mn)./e;

  % Compute real y-coordinates (obs. reversed y-axis direction)
  y=-S*(1-GAP)+repmat(ycoord,1,N)+(1-GAP)/2;          
  h_=plot(x',y');
  axis normal; 
  axis on; 
  set(gca,'ytick',1:M,'ydir','reverse'); 
  axis([0 N+1 0.5 M+.5]); grid off;
else
  h_=plot(S); 
  grid off;
end
set(h_,'color','k');

if nargout==1,
  h_line=h;
end
