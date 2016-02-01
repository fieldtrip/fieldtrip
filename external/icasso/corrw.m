function R=corrw(W,D)
%function R=corrw(W,D)
%
%PURPOSE
%
%To compute mutual linear correlation coefficients between M
%independent component estimates using the demixing matrix W and
%the dewhitening matrix D of the original data.  
%
% R=W*D*D'*W'; 
%
%INPUT
%
% W (field W in Icasso struct: demixing matrices)
% D (field dewhiteningMatrix in Icasso struct) the dewhitening
% matrix of the data 
%
%OUTPUT
%
% R (MxM matrix) of correlation coefficients
%
%SEE ALSO
% icassoSimilarity

% The normalization (rownorm) is a security measure: in some
% versions the FastICA has not normalized W properly if iteration
% stops prematurely

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

% ver 1.2 100105

B=rownorm(W*D);
R=B*B';

function X=rownorm(X)
% normalize rows to unit length.
s=abs(sqrt(sum(X.^2,2)));

if any(s==0),
  warning('Contains zero vectors: can''t normalize them!');
end

X=X./repmat(s,1,size(X,2));

