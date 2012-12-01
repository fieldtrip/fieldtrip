function X=parallelize(X,c)
%function X=parallelize(X,c)
%
%PURPOSE
%
%To multiply the columns of X by 1 or -1 as necessary so that the
%dot-products of column vectors of X are all positive (or zero)
%with respect to a reference vector c.     
%
%INPUT
%
% X (matrix) of size dxN of N d dimensional vectors
% c (vector) of size dx1, the reference vector 
%
%OUTPUT 
%
% X (matrix) of same size where the direction of some vectors has
% been reversed if necessary. 

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
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
%02111-1307, USA.

% ver 1.21 johan 020305

sgn=sign(c'*X);
sgn(sgn==0)=1;

X=X.*repmat(sgn,size(X,1),1);